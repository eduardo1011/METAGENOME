#!/usr/bin/env python
# coding: utf-8





import subprocess, warnings, os, shutil, re, pandas
warnings.filterwarnings("ignore")
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame
from functools import reduce
import math
import datetime

from ipywidgets import Button, GridBox, Layout, ButtonStyle, Box, HBox, VBox
import ipywidgets as widgets
# https://fontawesome.com/v4.7.0/icons/
from IPython.display import clear_output, display
import json
from scipy import stats
from collections import Counter, OrderedDict
from matplotlib.patches import Rectangle, Circle, Ellipse
import seaborn as sns
from matplotlib.colors import to_rgba_array, to_rgba, to_hex, to_rgb





# formatos de imagen= eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff





class HyperlinkManager:

    def __init__(self, text):

        self.text = text

        self.text.tag_config("hyper", foreground="blue", underline=1)

        self.text.tag_bind("hyper", "<Enter>", self._enter)
        self.text.tag_bind("hyper", "<Leave>", self._leave)
        self.text.tag_bind("hyper", "<Button-1>", self._click)

        self.reset()

    def reset(self):
        self.links = {}

    def add(self, action):
        # add an action to the manager.  returns tags to use in
        # associated text widget
        tag = "hyper-%d" % len(self.links)
        self.links[tag] = action
        return "hyper", tag

    def _enter(self, event):
        self.text.config(cursor="hand2")

    def _leave(self, event):
        self.text.config(cursor="")

    def _click(self, event):
        for tag in self.text.tag_names(CURRENT):
            if tag[:6] == "hyper-":
                self.links[tag]()
                return





os.makedirs('plots_asv',exist_ok=True)

os.makedirs('plots_asv/rarefaction',exist_ok=True)
os.makedirs('plots_asv/indices',exist_ok=True)
os.makedirs('plots_asv/richness',exist_ok=True)
os.makedirs('plots_asv/taxonomy',exist_ok=True)
os.makedirs('plots_asv/beta_diversity',exist_ok=True)











import matplotlib as mpl
from matplotlib import cm
import matplotlib
tab20 = [matplotlib.colors.to_hex(i) for i in cm.tab20(np.arange(20)/20.)]
tab10 = [matplotlib.colors.to_hex(i) for i in cm.tab10(np.arange(10)/10.)]
tab20b = [matplotlib.colors.to_hex(i) for i in cm.tab20b(np.arange(20)/20.)]
tab20c = [matplotlib.colors.to_hex(i) for i in cm.tab20c(np.arange(20)/20.)]
Set3 = [matplotlib.colors.to_hex(i) for i in cm.Set3(np.arange(12)/12.)]
Set2 = [matplotlib.colors.to_hex(i) for i in cm.Set2(np.arange(8)/8.)]
Set1 = [matplotlib.colors.to_hex(i) for i in cm.Set1(np.arange(9)/9.)]
Pastel2 = [matplotlib.colors.to_hex(i) for i in cm.Pastel2(np.arange(8)/8.)]
Pastel1 = [matplotlib.colors.to_hex(i) for i in cm.Pastel1(np.arange(9)/9.)]
Dark2 = [matplotlib.colors.to_hex(i) for i in cm.Dark2(np.arange(8)/8.)]
Paired = [matplotlib.colors.to_hex(i) for i in cm.Paired(np.arange(12)/12.)]
Accent = [matplotlib.colors.to_hex(i) for i in cm.Accent(np.arange(8)/8.)]
Spectral = [matplotlib.colors.to_hex(i) for i in cm.Spectral(np.arange(11)/11.)]
#Colormap1 = set(tab10 + Dark2 + Set1 + Paired + Set2 + Accent + tab20b  + Set3 + Pastel2 + tab20c+ Pastel1  + Spectral)

qualitative_colors = {'Pastel1':9, 'Pastel1_r':9,
                      'Pastel2':8, 'Pastel2_r':8,
                      'Paired':12, 'Paired_r':12,
                      'Accent':8, 'Accent_r':8,
                      'Dark2':8, 'Dark2_r':8,
                      'Set1':9, 'Set1_r':9,
                      'Set2':8, 'Set2_r':8,
                      'Set3':12, 'Set3_r':12,
                      'tab10':10, 'tab10_r':10,
                      'tab20':20, 'tab20_r':20,
                      'tab20b':20, 'tab20b_r':20,
                      'tab20c':20, 'tab20c_r':20}


qualitative_colors = {'Pastel1':9,
                      'Pastel2':8,
                      'Paired':12,
                      'Accent':8,
                      'Dark2':8,
                      'Set1':9,
                      'Set2':8,
                      'Set3':12,
                      'tab10':10,
                      'tab20':20,
                      'tab20b':20,
                      'tab20c':20}





column_scores = ['K_Score', 'P_Score', 'C_Score', 'O_Score', 'F_Score', 'G_Score', 'S_Score']
category_names = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']





metadata = pd.read_csv('METADATA.txt', sep = '\t')
metadata = metadata.astype('str')
metadata = metadata[metadata['Name Sample'].str.contains('ITS') == False]
metadata = metadata.sort_values(by =['Name Sample'],ascending=True).reset_index(drop=True)











SraRunTable = pd.read_csv('SraRunTable.txt', sep = ',')[['Sample Name', 'BioSample', 'Experiment']]
sample_BioSample = dict(zip(SraRunTable['Sample Name'], SraRunTable['BioSample']))











VARIABLE_KIT = 'Genomic DNA kit'





DNPowerSoil = metadata[metadata[VARIABLE_KIT] == 'DNPowerSoil']['Name Sample'].tolist()
DNMicrobial = metadata[metadata[VARIABLE_KIT] == 'DNMicrobial']['Name Sample'].tolist()
name_sample = metadata['Name Sample'].tolist()





variables = ['Coffee Variety', 'Cultivation', 'Postharvest Processing', VARIABLE_KIT,
             'Drying Time (Days)', 'ug OTA/kg', 'Location']
variables = list(reversed(list(dict(OrderedDict(Counter({i : len(i) for i in variables}).most_common())).keys())))





correspondencia_sam_vars = {}
for i, row in metadata[['Name Sample'] + variables].iterrows():
    correspondencia_sam_vars[row['Name Sample']] = {'Genomic DNA kit':row[VARIABLE_KIT], 'Location':row['Location'], 'ug OTA/kg':row['ug OTA/kg'],
                                                    'Cultivation':row['Cultivation'], 'Coffee Variety':row['Coffee Variety'],
                                                    'Drying Time (Days)':row['Drying Time (Days)'], 'Postharvest Processing':row['Postharvest Processing']}











name_code = dict(zip(metadata['Name Sample'], metadata['Code']))
name_code2 = dict(zip(metadata['Name Sample'], metadata['Code2']))
code_name = dict(zip(metadata['Code'], metadata['Name Sample']))
code2_name = dict(zip(metadata['Code2'], metadata['Name Sample']))





KITS = {'Both kits' : name_sample}





for i in metadata[VARIABLE_KIT].unique():
    lista = metadata[metadata[VARIABLE_KIT] == i]['Name Sample'].tolist()
    KITS.update({i : lista})











"""
# ultra FAST
Créditos a: [Quick and dirty way to calculate large binomial coefficients in Python]
(https://grocid.net/2012/07/02/quick-and-dirty-way-to-calculate-large-binomial-coefficients-in-python/)
"""
def log_fac2(m,n):
    from math import log
    return n*log(n)-m*log(m)
 
def log_binomial2(n,k):
    if k > (n-k):
        return (log_fac2(n-k+1,n)-log_fac2(2,k))
    else:
        return (log_fac2(k+1,n)-log_fac2(2,n-k))











marker_dict = ['o', 'v', '^', '<', '>', 's', 'p', 'H', '*', 'D']
width_line = []
n = 0.5
for i in range(20):
    width_line.append(n)
    n += 0.5





ASV_counts = pd.read_csv('clustering/ASV_counts.txt', sep = '\t')


#### base de datos

# Selections
Columnas = [0, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23, 25]
Nombres = ['Entry', 'Kingdom', 'K_Score', 'Phylum', 'P_Score', 'Class', 'C_Score',
           'Order', 'O_Score', 'Family', 'F_Score', 'Genus', 'G_Score', 'Species', 'S_Score']
column_scores = ['K_Score', 'P_Score', 'C_Score', 'O_Score', 'F_Score', 'G_Score', 'S_Score']
category_names = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
lin_scores = dict(zip(category_names, column_scores))
scores_lin = dict(zip(column_scores, category_names))
asv_classif, asv = ['ASV_vs_NCBI_classif', 'ASV_vs_RDP_classif', 'ASV_vs_SILVA_classif'], 'ASV'

NCBI = {}
RDP = {}
SILVA = {}
for _clas in asv_classif:
    df = pd.read_csv('tablas/'+_clas+'.txt', sep = '\t', header = None)
    df = df[Columnas]
    df.columns = Nombres
    df = df[df.Species.str.contains('virus|_strain') == False]
    if 'NCBI' == _clas.split('_')[2]:
        for e, i in enumerate(df['Entry']):
            NCBI[i] = list(df.iloc[e].values[1:])
    if 'RDP' == _clas.split('_')[2]:
        for e, i in enumerate(df['Entry']):
            RDP[i] = list(df.iloc[e].values[1:])
    if 'SILVA' == _clas.split('_')[2]:
        for e, i in enumerate(df['Entry']):
            SILVA[i] = list(df.iloc[e].values[1:])
colnames = df.columns[1:].tolist()




help0_button = widgets.Button(description="Help", icon = 'fa-question-circle', layout=Layout(width='65px'))
help0_button.style.button_color = 'pink'
help0_button.style.font_weight = 'bold'
help0_button_output = widgets.Output()

import tkinter as tk
from tkinter import *
    
def button_clicked(b):
    with help0_button_output:
        clear_output(True)
        
        import webbrowser
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)

        newwin0 = Tk()
        newwin0.title("Information about analysis")

        newwin0.geometry("515x200")
        newwin0.configure(background='white')
        newwin0.resizable(0, 0)
        newwin0.attributes('-topmost', 1)

        c00 = Label(newwin0, text="    ", bg = 'white')
        c00.grid(column=0, row=0)

        text2 = Text(newwin0, height=40, width=80, bg = 'white', bd = 0, font=('Arial', 7))
        scroll = Scrollbar(newwin0, command=text2.yview)
        text2.configure(yscrollcommand=scroll.set)

        text2.tag_configure('big', font=('Arial', 8, 'bold'))
        # titulos

        # contenidos
        text2.tag_configure('color', foreground='red', font=('Arial', 7, 'bold'))
        text2.tag_configure('color2', foreground='darkorange', font=('Arial', 7, 'bold'))

        text2.tag_configure('text', font=('Arial', 7))
        text2.tag_configure('text2', font=('Arial',7, 'bold'))

        hyperlink = HyperlinkManager(text2)

        #------------------------------------------------------------
        text2.insert(END,'\nScores assigned based on the RDP ranking algorithm.\n\n', 'big')
        
        def click1():
            webbrowser.open_new(r"https://raw.githubusercontent.com/eduardo1011/Bioinformatica2019/master/RDP_algorithm.JPG")
        text2.insert(END, "Look at the following example.\n\n", hyperlink.add(click1))
        
        text2.insert(END, "\n\nFor more information on the RDP algorithm visit:\n", 'text')
        def click1():
            webbrowser.open_new(r"http://rdp.cme.msu.edu/")
        text2.insert(END, "http://rdp.cme.msu.edu/\n\n", hyperlink.add(click1))
        
        def click1():
            webbrowser.open_new(r"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950982/")
        text2.insert(END, "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950982/\n\n", hyperlink.add(click1))
        
        text2.grid(column=1, row=0)
        scroll.grid(column=1, row=0, sticky = 'ESN')
        text2.configure(state=DISABLED)
        newwin0.mainloop()
        
         
help0_button.on_click(button_clicked)
ayuda0 = HBox([help0_button, help0_button_output])



blanca = widgets.Button(layout=Layout(width='10px', height='25px'), disabled=True)
blanca.style.button_color = 'white'

progress = widgets.HTML()

tipo_bd = widgets.ToggleButtons(options= ['NCBI', 'RDP', 'SILVA', 'MERGE ALL DBs'], value = 'NCBI', button_style = '') # este es solo para 16S
tipo_bd.style.button_width = '150px'
tipo_bd.style.font_weight = 'bold'

tipo_bd_box = Box(children=[VBox([tipo_bd])], layout=Layout(border='1px solid silver', width='160px', height='120px'))

umbral_tax_genus = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.5,description='Genus >=:',disabled=False,
                                continuous_update=False,orientation='horizontal',readout=True,
                                   layout=Layout(width='300px', height='25px'))
umbral_tax_specie = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.1,description='Species >=:',disabled=False,
                                continuous_update=False,orientation='horizontal',readout=True,
                                   layout=Layout(width='300px', height='25px'))
limite_reads = widgets.SelectionSlider(options=range(0, 201),value=2,description='Reads >=:',disabled=False,
                                continuous_update=False,orientation='horizontal',readout=True,
                                   layout=Layout(width='300px', height='25px'))


confirmacion = widgets.ToggleButtons(options=['False', 'True'])
confirmacion.style.button_width = '50px'



def box1(tipo_bd):
    
        
    if tipo_bd == 'MERGE ALL DBs':
        progress.value = '<font color = black> <i class="fa fa-spinner fa-pulse fa-3x fa-fw"></i> </font>'
        #---------------
        result = []
        for e, entry in enumerate(set(list(NCBI.keys()) + list(RDP.keys()) + list(SILVA.keys()))):
            record = []
            GGG = []
            for db in [NCBI, RDP, SILVA]:
                if entry in list(db.keys()):
                    valores = [float(db[entry][e+1]) for e, i in enumerate(db[entry]) if e%2 == 0]
                    record.append([i for e, i in enumerate(db[entry]) if e%2 == 0] + valores + [np.mean(valores[4:])])
            record.sort(key = lambda x: x[14], reverse = True)
            result.append([entry] + record[0])
        NCBI_RDP_SILVA = DataFrame(result, columns = ['Entry'] + category_names + column_scores + ['Mean'])
        #-----------------
        
        Unification  = []
        for i in NCBI_RDP_SILVA.Species.drop_duplicates():
            """
            toma la media mas alta para especies duplicadas con diferente linaje taxonomico
            """
            df = NCBI_RDP_SILVA[NCBI_RDP_SILVA.Species == i].sort_values(by ='Mean',ascending=False).reset_index(drop=True).sort_values(by ='Mean',ascending=False).reset_index(drop=True)
            Unification.append(df[category_names].values.tolist()[0])
            
        ASV_Full_Taxonomy = NCBI_RDP_SILVA[['Entry', 'Species'] + column_scores + ['Mean']].merge(DataFrame(Unification, columns = category_names), on = 'Species', how = 'left')
        ASV_Full_Taxonomy = ASV_Full_Taxonomy.merge(ASV_counts, on = 'Entry', how = 'left')
        
        ASV_Full_Taxonomy.to_csv('tablas/ASV_Full_Taxonomy.txt', sep = '\t', index = None)
        
        clear_output(True)
        progress.value = '<font color = black> <i class="fa fa-spinner fa-3x fa-fw"></i> </font>'
        
    if tipo_bd == 'NCBI': 
        progress.value = '<font color = black> <i class="fa fa-spinner fa-pulse fa-3x fa-fw"></i> </font>'
        #---------------
        result = []
        for e, entry in enumerate(set(list(NCBI.keys()))):
            record = []
            GGG = []
            for db in [NCBI]:
                if entry in list(db.keys()):
                    valores = [float(db[entry][e+1]) for e, i in enumerate(db[entry]) if e%2 == 0]
                    record.append([i for e, i in enumerate(db[entry]) if e%2 == 0] + valores + [np.mean(valores[4:])])
            record.sort(key = lambda x: x[14], reverse = True)
            result.append([entry] + record[0])
        NCBI_DB = DataFrame(result, columns = ['Entry'] + category_names + column_scores + ['Mean'])
        #---------------
        
        Unification  = []
        for i in NCBI_DB.Species.drop_duplicates():
            """
            toma la media mas alta para especies duplicadas con diferente linaje taxonomico
            """
            df = NCBI_DB[NCBI_DB.Species == i].sort_values(by ='Mean',ascending=False).reset_index(drop=True).sort_values(by ='Mean',ascending=False).reset_index(drop=True)
            Unification.append(df[category_names].values.tolist()[0])
            
        ASV_Full_Taxonomy = NCBI_DB[['Entry', 'Species'] + column_scores + ['Mean']].merge(DataFrame(Unification, columns = category_names), on = 'Species', how = 'left')
        ASV_Full_Taxonomy = ASV_Full_Taxonomy.merge(ASV_counts, on = 'Entry', how = 'left')
        ASV_Full_Taxonomy.to_csv('tablas/ASV_Full_Taxonomy.txt', sep = '\t', index = None)
        
        clear_output(True)
        progress.value = '<font color = black> <i class="fa fa-spinner fa-3x fa-fw"></i> </font>'
        
    if tipo_bd == 'RDP':
        progress.value = '<font color = black> <i class="fa fa-spinner fa-pulse fa-3x fa-fw"></i> </font>'
        #---------------
        result = []
        for e, entry in enumerate(set(list(RDP.keys()))):
            record = []
            GGG = []
            for db in [RDP]:
                if entry in list(db.keys()):
                    valores = [float(db[entry][e+1]) for e, i in enumerate(db[entry]) if e%2 == 0]
                    record.append([i for e, i in enumerate(db[entry]) if e%2 == 0] + valores + [np.mean(valores[4:])])
            record.sort(key = lambda x: x[14], reverse = True)
            result.append([entry] + record[0])
        RDP_DB = DataFrame(result, columns = ['Entry'] + category_names + column_scores + ['Mean'])
        #---------------
        
        Unification  = []
        for i in RDP_DB.Species.drop_duplicates():
            """
            toma la media mas alta para especies duplicadas con diferente linaje taxonomico
            """
            df = RDP_DB[RDP_DB.Species == i].sort_values(by ='Mean',ascending=False).reset_index(drop=True).sort_values(by ='Mean',ascending=False).reset_index(drop=True)
            Unification.append(df[category_names].values.tolist()[0])
            
        ASV_Full_Taxonomy = RDP_DB[['Entry', 'Species'] + column_scores + ['Mean']].merge(DataFrame(Unification, columns = category_names), on = 'Species', how = 'left')
        ASV_Full_Taxonomy = ASV_Full_Taxonomy.merge(ASV_counts, on = 'Entry', how = 'left')
        ASV_Full_Taxonomy.to_csv('tablas/ASV_Full_Taxonomy.txt', sep = '\t', index = None)
        
        clear_output(True)
        progress.value = '<font color = black> <i class="fa fa-spinner fa-3x fa-fw"></i> </font>'
        
        
    if tipo_bd == 'SILVA':
        progress.value = '<font color = black> <i class="fa fa-spinner fa-pulse fa-3x fa-fw"></i> </font>'
        #---------------
        result = []
        for e, entry in enumerate(set(list(SILVA.keys()))):
            record = []
            GGG = []
            for db in [SILVA]:
                if entry in list(db.keys()):
                    valores = [float(db[entry][e+1]) for e, i in enumerate(db[entry]) if e%2 == 0]
                    record.append([i for e, i in enumerate(db[entry]) if e%2 == 0] + valores + [np.mean(valores[4:])])
            record.sort(key = lambda x: x[14], reverse = True)
            result.append([entry] + record[0])
        SILVA_DB = DataFrame(result, columns = ['Entry'] + category_names + column_scores + ['Mean'])
        #---------------
        
        Unification  = []
        for i in SILVA_DB.Species.drop_duplicates():
            """
            toma la media mas alta para especies duplicadas con diferente linaje taxonomico
            """
            df = SILVA_DB[SILVA_DB.Species == i].sort_values(by ='Mean',ascending=False).reset_index(drop=True).sort_values(by ='Mean',ascending=False).reset_index(drop=True)
            Unification.append(df[category_names].values.tolist()[0])
            
        ASV_Full_Taxonomy = SILVA_DB[['Entry', 'Species'] + column_scores + ['Mean']].merge(DataFrame(Unification, columns = category_names), on = 'Species', how = 'left')
        ASV_Full_Taxonomy = ASV_Full_Taxonomy.merge(ASV_counts, on = 'Entry', how = 'left')
        ASV_Full_Taxonomy.to_csv('tablas/ASV_Full_Taxonomy.txt', sep = '\t', index = None)
        
        clear_output(True)
        progress.value = '<font color = black> <i class="fa fa-spinner fa-3x fa-fw"></i> </font>'
        
        
        
    ############################
    
    conteo_inicial = VBox([widgets.HTML('<font color = grey> <b style="font-size:0.6vw">ASVs : '+str(len(set(ASV_Full_Taxonomy.Entry.unique())))+'</b>'),
                           widgets.HTML('<font color = grey> <b style="font-size:0.6vw">Genus : '+str(len(set(ASV_Full_Taxonomy.Genus.unique())))+'</b>'),
                           widgets.HTML('<font color = grey> <b style="font-size:0.6vw">Species : '+str(len(set(ASV_Full_Taxonomy.Species.unique())))+'</b>')])
    conteo_inicial_box = Box(children=[conteo_inicial], layout=Layout(display='flex', flex_flow='column', align_items='stretch', border='1px solid grey',
                                                                      width='120px', height=str(int(len(conteo_inicial.children) * 34))+'px'))
    
    
    
    
    def box2(umbral_tax_genus, umbral_tax_specie, limite_reads):    
    
    
        Full_Taxonomy = ASV_Full_Taxonomy[ASV_Full_Taxonomy.G_Score >= umbral_tax_genus] # threshold of 0.8 at the genus level
        Full_Taxonomy = Full_Taxonomy[Full_Taxonomy.S_Score >= umbral_tax_specie] # threshold of 0.5 at the species level

        COUNTS_THRESDHOLD = limite_reads
        filtrados = {i : Full_Taxonomy[['Entry', i]].values[Full_Taxonomy[['Entry', i]].values[:, 1] >= COUNTS_THRESDHOLD] for i in name_sample}
        frames = [DataFrame(filtrados[i].tolist(), columns = ['Entry', i]) for i in filtrados]
        NCBI_RDP_SILVA_SUMMARY = reduce(lambda  left,right: pd.merge(left, right, on = ['Entry'], how = 'outer'), frames).fillna(float(0))
        especies = DataFrame(NCBI_RDP_SILVA_SUMMARY.Entry.tolist(), columns = ['Entry']).merge(ASV_Full_Taxonomy[['Entry', 'Species']], on = 'Entry', how = 'left').Species.tolist()
        generos = DataFrame(NCBI_RDP_SILVA_SUMMARY.Entry.tolist(), columns = ['Entry']).merge(ASV_Full_Taxonomy[['Entry', 'Genus']], on = 'Entry', how = 'left').Genus.tolist()
        NCBI_RDP_SILVA_SUMMARY.insert(loc = 1, column='Species', value = especies)
        NCBI_RDP_SILVA_SUMMARY.insert(loc = 2, column='Genus', value = generos)

        #print('  Result')
        #print('   ASVs:', len(set(NCBI_RDP_SILVA_SUMMARY.Entry.unique())))
        #print('Species:', len(set(NCBI_RDP_SILVA_SUMMARY.Species.unique())))
        #print('  Genus:', len(set(NCBI_RDP_SILVA_SUMMARY.Genus.unique())))


        procesar_button = widgets.Button(description="PROCESS", icon = 'fa-play', layout=Layout(width='120px'))
        procesar_button.style.button_color = 'lime' #'deepskyblue'
        procesar_button.style.font_weight = 'bold'
        procesar_output = widgets.Output()


        resultado = VBox([widgets.HTML('<b style="font-size:0.6vw">ASVs : '+str(len(set(NCBI_RDP_SILVA_SUMMARY.Entry.unique())))+'</b>'),
                           widgets.HTML('<b style="font-size:0.6vw">Genus : '+str(len(set(NCBI_RDP_SILVA_SUMMARY.Genus.unique())))+'</b>'),
                          widgets.HTML('<b style="font-size:0.6vw">Species : '+str(len(set(NCBI_RDP_SILVA_SUMMARY.Species.unique())))+'</b>')
                           ])
        resultado_box = Box(children=[resultado], layout=Layout(display='flex', flex_flow='column', align_items='stretch', border='3px solid limegreen',
                                                                width='120px', height=str(int(len(conteo_inicial.children) * 34))+'px'))



        #pro = HBox([VBox([blanca, procesar_button]),
        #           VBox([blanca, widgets.HTML('<font color = black> <i class="fa fa-play fa-2x" aria-hidden="true"></i>')])])

        #display(HBox([pro, blanca, resultado_box]))

        display(HBox([conteo_inicial_box, blanca, resultado_box, blanca, ayuda0, VBox([procesar_button, procesar_output])]))



        def button_clicked(b):
            with procesar_output:
                clear_output(True)



                NCBI_RDP_SILVA_SUMMARY.to_csv('tablas/ASVs_NCBI_RDP_SILVA_SUMMARY.txt',index = None, sep = '\t')
                
                Full_Taxonomy.to_csv('tablas/ASV_Full_Taxonomy_Threshold.txt', sep = '\t', index = None)
                datasets = []
                for m in name_sample:
                    df2 = Full_Taxonomy[category_names+[m]][Full_Taxonomy[category_names+[m]][m] >= COUNTS_THRESDHOLD].sort_values(by =m,ascending=False).reset_index(drop=True)
                    df2 = df2[[m]+category_names]
                    cuentas = df2[m].sum()
                    df2['Ratio'] = (df2[m] / cuentas) * 100
                    df2 = df2.rename(columns={m: 'Counts'})
                    df2.insert(loc = 0, column='Sample', value=[m]*len(df2))
                    datasets.append(df2)
                DataSets = pd.concat(datasets)
                """
                Dataframe generado con datos para mostrar con el entorno interactivo de Krona
                """
                DataSets[['Sample', 'Counts', 'Ratio'] + category_names].to_csv('tablas/ASV_Full_Taxonomy_for_KRONA.txt',index = None, sep = '\t', header = None)


                #display(resultado_box)

                #rarefaction = NCBI_RDP_SILVA_SUMMARY[['Entry'] + name_sample]
                rarefaction_collapse_specie = pd.pivot_table(NCBI_RDP_SILVA_SUMMARY[['Species'] + name_sample], values = name_sample, index = ['Species'], aggfunc = sum).reset_index()
                #rarefaction_collapse_specie.to_csv('tablas/ASVs_rarefaction_collapse_specie.txt',index = None, sep = '\t')


                #data_for_rarefaction = {i : [{'No_species': len(np.array(rarefaction[i])[np.array(rarefaction[i]) > 0].tolist())},
                #                             {'Counts': np.array(rarefaction[i])[np.array(rarefaction[i]) > 0].tolist()}] for i in name_sample}
                data_for_rarefaction = {i : [{'No_species': len(np.array(rarefaction_collapse_specie[i])[np.array(rarefaction_collapse_specie[i]) > 0].tolist())},
                                             {'Counts': np.array(rarefaction_collapse_specie[i])[np.array(rarefaction_collapse_specie[i]) > 0].tolist()}] for i in name_sample}

                with open('tablas/ASVs_data_for_rarefaction.json', 'w') as fp:
                    json.dump(data_for_rarefaction, fp)

                procesando = widgets.IntProgress(value=0, min=0, max=10, description='Processing:', bar_style='', style={'bar_color': 'black'}, orientation='horizontal', layout=Layout(width='300px', height='25px'))

                import time
                display(procesando)
                for k in range(10):
                    time.sleep(0.05)
                    procesando.value = k+1

                clear_output(True)
                #----------------------------------------------
        procesar_button.on_click(button_clicked)

    OUT2 = widgets.interactive_output(box2, {'umbral_tax_genus':umbral_tax_genus, 'umbral_tax_specie':umbral_tax_specie, 'limite_reads':limite_reads})
    
    display(OUT2)

OUT1 = widgets.interactive_output(box1, {'tipo_bd':tipo_bd})


threshold = HBox([HBox([VBox([widgets.HTML('<font color = grey><b style="font-size:0.8vw">SELECT A DATABASE: </b>'), progress]), tipo_bd_box]),
                  blanca, widgets.HTML('<font color = black> <i class="fa fa-filter fa-2x" aria-hidden="true"></i> <h style="font-size:0.7vw">Taxonomic threshold:</h>'),
                  VBox([umbral_tax_genus, umbral_tax_specie, limite_reads]), OUT1])





# # Rarefaction




help3_button = widgets.Button(description="Help", icon = 'fa-question-circle', layout=Layout(width='65px'))
help3_button.style.button_color = 'pink'
help3_button.style.font_weight = 'bold'
help3_button_output = widgets.Output()

    
def button_clicked(b):
    with help3_button_output:
        clear_output(True)
        
        import webbrowser
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)

        newwin0 = Tk()
        newwin0.title("Information about Colors")

        newwin0.geometry("515x150")
        newwin0.configure(background='white')
        newwin0.resizable(0, 0)
        newwin0.attributes('-topmost', 1)

        c00 = Label(newwin0, text="    ", bg = 'white')
        c00.grid(column=0, row=0)

        text2 = Text(newwin0, height=40, width=80, bg = 'white', bd = 0, font=('Arial', 7))
        scroll = Scrollbar(newwin0, command=text2.yview)
        text2.configure(yscrollcommand=scroll.set)

        text2.tag_configure('big', font=('Arial', 8, 'bold'))
        # titulos

        # contenidos
        text2.tag_configure('color', foreground='red', font=('Arial', 7, 'bold'))
        text2.tag_configure('color2', foreground='darkorange', font=('Arial', 7, 'bold'))

        text2.tag_configure('text', font=('Arial', 7))
        text2.tag_configure('text2', font=('Arial',7, 'bold'))

        hyperlink = HyperlinkManager(text2)

        #------------------------------------------------------------
        text2.insert(END,'\nColors\n', 'big')
        
        text2.insert(END, '\nYou can select more than one code from the list to assign a color to each sample.\n', 'text')
        
        text2.insert(END, "\nSee Qualitative Colormaps in the link:\n", 'text')
        
        def click1():
            webbrowser.open_new(r"https://matplotlib.org/stable/tutorials/colors/colormaps.html")
        text2.insert(END, "https://matplotlib.org/stable/tutorials/colors/colormaps.html\n\n", hyperlink.add(click1))
        
        
        text2.grid(column=1, row=0)
        scroll.grid(column=1, row=0, sticky = 'ESN')
        text2.configure(state=DISABLED)
        newwin0.mainloop()
        
         
help3_button.on_click(button_clicked)
ayuda3 = HBox([help3_button, help3_button_output])











"""
Rarefaction analysis
"""
from matplotlib.colors import ListedColormap # funcion para crear un objeto <matplotlib.colors.ListedColormap> a partir de una lista de colores personalizados


from tkinter import ttk

progreso1 = widgets.HTML('<font color = black> <i class="fa fa-spinner fa-2x fa-fw"></i> </font>')
progreso2 = widgets.HTML()
estatico = widgets.HTML('<font color = limegreen> <b style="font-size:0.5vw">Processed samples : </b>')
estatico = HBox([progreso1, estatico, progreso2])
estatico_box = Box(children=[estatico], layout=Layout(border='1px solid limegreen', width='350px', height='32px'))

RARE_button = widgets.Button(description="PROCESS AND VISUALIZE", icon = 'fa-eye', layout=Layout(width='590px'))
RARE_button.style.button_color = 'gainsboro' #'deepskyblue'
RARE_button.style.font_weight = 'bold'
RARE_output = widgets.Output()

def button_clicked(b):
    import time
    import ctypes
    
    with open('tablas/ASVs_data_for_rarefaction.json', 'r') as fp:
        data_for_rarefaction = json.load(fp)
    
    with RARE_output:
        clear_output(True)
        progreso1.value = '<font color = black> <i class="fa fa-spinner fa-pulse fa-2x fa-fw"></i> </font>'
        #progreso = widgets.HTML()
        #progress = HBox([widgets.HTML('<font color = black> <i class="fa fa-spinner fa-pulse fa-2x fa-fw"></i> </font>'), progreso])
        #display(Box(children=[progress], layout=Layout(border='1px solid limegreen', width='585px', height='32px')))
        
        
        RAREFACTION = {}
        for G in data_for_rarefaction:
            
            #progreso.value = '<font color = limegreen> <b style="font-size:0.5vw">Processed samples : </b> <font color = black> <b style="font-size:0.5vw"> '+G+'</b>'
            progreso2.value = '</b> <font color = black> <b style="font-size:0.5vw"> '+G+'</b>'
            
            Counts = data_for_rarefaction[G][1]['Counts']
            SUMA = np.sum(Counts)
            #print(G, len(Counts), SUMA)
            sample = np.arange(1, int(SUMA), 100)
            xdata = []
            ydata = []
            for s in sample:
                ldiv = log_binomial2(SUMA,s)
                U = []
                for h in [int(i) for i in list(SUMA - np.array(Counts, dtype = 'int'))]:
                    if h < s:
                        U.append(0)
                    else:
                        U.append(log_binomial2(int(h), s))
                p1 = np.exp(np.array(U) - ldiv)
                out = np.sum(1 - p1)
                xdata.append(s)
                ydata.append(out)
            RAREFACTION[G] = {'xdata':np.array(xdata), 'ydata':np.array(ydata)}
            time.sleep(0.35)

        val_min_x = min(set([j for i in RAREFACTION for j in RAREFACTION[i]['xdata']]))
        val_max_x = max(set([j for i in RAREFACTION for j in RAREFACTION[i]['xdata']]))
        val_min_y = min(set([j for i in RAREFACTION for j in RAREFACTION[i]['ydata']]))
        val_max_y = max(set([j for i in RAREFACTION for j in RAREFACTION[i]['ydata']]))
        
        
        clear_output(True)
        progreso1.value = '<font color = black> <i class="fa fa-spinner fa-2x fa-fw"></i> </font>'
        progreso2.value = '</b> <font color = black> <b style="font-size:0.5vw"> '+G+'</b>'

        #progress = HBox([widgets.HTML('<font color = black> <i class="fa fa-spinner fa-2x fa-fw"></i> </font>'), progreso])
        #display(Box(children=[progress], layout=Layout(border='1px solid limegreen', width='585px', height='32px')))
        
        
        mpl.rcParams.update(mpl.rcParamsDefault)

        fig1 = plt.figure(figsize=(6, 4))

        AX1 = fig1.add_axes([0, 0, 1, 1])

        plt.gca().tick_params(which='major', width = 2, length=4, color='gainsboro')
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_color('gainsboro')
        plt.gca().spines['bottom'].set_color('gainsboro')
        plt.gca().spines['right'].set_color(None)
        plt.gca().spines['top'].set_color(None)

        plt.gca().set_xlabel('Contigs', fontsize=12, fontname='Open Sans', weight = 'bold')
        plt.gca().set_ylabel('Species', fontsize=12, fontname='Open Sans', weight = 'bold')
        plt.close()
        
        
        fig_size = []
        n = 3
        for i in range(16):
            fig_size.append(round(n, 2))
            n += 0.2
    
        tam_plot1 = widgets.SelectionSlider(options=fig_size,value=4,disabled=False,
                                              description = 'Chart size:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        Multiple_colorS = widgets.SelectMultiple(options=sorted(list(qualitative_colors.keys())), value=sorted(list(qualitative_colors.keys()))[0:5], disabled=False,
                       layout=Layout(width='110px', height='250px'))
        ## lineas
        width_linea = widgets.SelectionSlider(options=width_line,value=1,disabled=False,
                                              description = 'Line width:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                             layout=Layout(width='400px', height='25px'))
        linea_color = widgets.ColorPicker(concise=True, value='black', disabled=False, layout = Layout(width='35px', height='25px'))
        
        alfa_linea = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=1,disabled=False,
                                              description = 'Line alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        # puntos

        size_point = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.5,disabled=False,
                                              description = 'Marker size:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        alfa_point = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=1,disabled=False,
                                              description = 'Marker alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))

        #markers_point = widgets.SelectionSlider(options= marker_dict,value='o',disabled=False,
        #                                      description = 'Marker:', continuous_update=False,orientation='horizontal',readout=True,
        #                                    layout=Layout(width='400px', height='25px'))
        markers_point = widgets.ToggleButtons(options=['⏺', '⏹'], value = '⏺')
        markers_point.style.button_width = '31px'
        ### leyenda
        size_title_legend = widgets.SelectionSlider(options=range(5, 21),value=10,disabled=False,
                                              description = 'Legend title:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        size_text_legend = widgets.SelectionSlider(options=range(5, 21),value=10,disabled=False,
                                              description = 'Legend text:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        family = sorted(['Liberation Serif','Microsoft Sans Serif','Open Sans','Times New Roman','3ds Light','Calibri','Comic Sans MS',
                  'Arial','Courier New','Microsoft Yi Baiti','Lucida Console'])
        family1 = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                   layout = Layout(width='290px', height='25px'))

        #num_cols = widgets.SelectionSlider(options=range(1, 6),value=1,disabled=False,
        #                                      description = 'Num cols:', continuous_update=False,orientation='horizontal',readout=True,
        #                                    layout=Layout(width='400px', height='25px'))

        num_cols = widgets.ToggleButtons(options=[1, 2, 3, 4], value = 2)
        num_cols.style.button_width = '25px'

        incluir_cuentas = widgets.ToggleButtons(options=['True', 'False'], value = 'False')
        incluir_cuentas.style.button_width = '55px'

        ### ejes
        label_X = widgets.SelectionSlider(options=range(5, 31),value=12,disabled=False,
                                              description = 'X label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_X = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'X tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))

        label_Y = widgets.SelectionSlider(options=range(5, 31),value=12,disabled=False,
                                              description = 'Y label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_Y = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'Y tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))

        family_axis = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                   layout = Layout(width='290px', height='25px'))


        cambiar_sample = widgets.ToggleButtons(options=['Code1', 'Code2', 'Code3'], value = 'Code1')
        cambiar_sample.style.button_width = '75px'
        
        ver_num_muestra = widgets.ToggleButtons(options=['True', 'False'], value = 'False')
        ver_num_muestra.style.button_width = '55px'
        
        ajustar_ejes = widgets.ToggleButtons(options=['True', 'False'], value = 'False')
        ajustar_ejes.style.button_width = '55px' 


        tipo_kits_1 = widgets.ToggleButtons(options= ['Both kits'] + metadata[VARIABLE_KIT].unique().tolist(), value = 'Both kits', button_style = 'primary')
        tipo_kits_1.style.button_width = '125px'
        tipo_kits_1.style.font_weight = 'bold'

        tipo_kit_1_box = Box(children=[VBox([tipo_kits_1])], layout=Layout(border='1px solid #1976d2', width='135px', height='95px'))
        
        ##
        filename_plot = widgets.Text(value='', placeholder='Chart name', description='', disabled=False, layout = Layout(width='270px', height='25px'))
        
        ancho_plot_save = str(71)
        png1 = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
        png1.style.button_color = 'gold'
        output1 = widgets.Output()
        def button_clicked1(b):
            with output1:
                clear_output(True)
                if filename_plot.value is '':
                    nombre_grafico = 'Rarefy_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                if filename_plot.value is not '':
                    nombre_grafico = filename_plot.value
                fig1.savefig('plots_asv/rarefaction/'+nombre_grafico+'.png', dpi = 900, bbox_inches= 'tight')
        png1.on_click(button_clicked1)
        #----
        #----
        jpeg1 = widgets.Button(description="JPEG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
        jpeg1.style.button_color = 'gold'
        output2 = widgets.Output()
        def button_clicked2(b):
            with output2:
                clear_output(True)
                if filename_plot.value is '':
                    nombre_grafico = 'Rarefy_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                if filename_plot.value is not '':
                    nombre_grafico = filename_plot.value
                fig1.savefig('plots_asv/rarefaction/'+nombre_grafico+'.jpeg', dpi = 900, bbox_inches= 'tight')
                
        jpeg1.on_click(button_clicked2)
        #----
        svg1 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
        svg1.style.button_color = 'gold'
        output3 = widgets.Output()
        def button_clicked3(b):
            with output3:
                clear_output(True)
                if filename_plot.value is '':
                    nombre_grafico = 'Rarefy_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                if filename_plot.value is not '':
                    nombre_grafico = filename_plot.value
                fig1.savefig('plots_asv/rarefaction/'+nombre_grafico+'.svg', dpi = 900, bbox_inches= 'tight')
                
        svg1.on_click(button_clicked3)
        #----
        pdf1 = widgets.Button(description="PDF", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
        pdf1.style.button_color = 'gold'
        output4 = widgets.Output()
        def button_clicked4(b):
            with output4:
                clear_output(True)
                if filename_plot.value is '':
                    nombre_grafico = 'Rarefy_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                if filename_plot.value is not '':
                    nombre_grafico = filename_plot.value
                fig1.savefig('plots_asv/rarefaction/'+nombre_grafico+'.pdf', dpi = 900, bbox_inches= 'tight')
                
        pdf1.on_click(button_clicked4)

        #----
        filename = widgets.Text(value='', placeholder='File name.txt', description='', disabled=False, layout = Layout(width='195px', height='25px'))
        
        params = widgets.Button(description="SAVE", icon = 'fa-save', layout=Layout(width='70px'))
        params.style.button_color = 'cyan'
        output5 = widgets.Output()
        def button_clicked5(b):
            with output5:
                clear_output(True)
                
                
                if filename.value is '':
                    nombre_archivo = 'Rarefy_chart_params_ASVs.txt'
                if filename.value is not '':
                    nombre_archivo = filename.value

                with open('plots_asv/rarefaction/'+nombre_archivo, 'w') as fq:
                    fq.write('#Saved parameters\n')
                    fq.write('Clustering:ASVs\n' )
                    """
                    diccionario para guardar en un archivo los parametros usados para volver a reproducirlos posteriormente
                    """
                    parametros = {'Kit:':tipo_kits_1.value,
                                  'Adjust axes:':ajustar_ejes.value,
                                  'Chart size:':tam_plot1.value,
                                  'Line width:':width_linea.value,
                                  'Line color':linea_color.value,
                                  'Line alpha:':alfa_linea.value,
                                  'Marker colors:':', '.join(list(Multiple_colorS.value)),
                                  'Marker size:':size_point.value,
                                  'Marker alpha:':alfa_point.value,
                                  'Marker:':{'⏺':'circle', '⏹':'square'}[markers_point.value],
                                  'Sample number:':ver_num_muestra.value,
                                  'X label:':label_X.value,
                                  'X tick labels:':ticklabels_X.value,
                                  'Y label:':label_Y.value,
                                  'Y tick labels:':ticklabels_Y.value,
                                  'Axis font:':family_axis.value,
                                  'Legend title:':size_title_legend.value,
                                  'Legend text:':size_text_legend.value,
                                  'Legend font:':family1.value,
                                  'Columns:':num_cols.value,
                                  'Contigs:':incluir_cuentas.value,
                                  'Rename:':cambiar_sample.value}
                    
                    for w in parametros:
                        fq.write(w+str(parametros[w])+'\n')
        
        params.on_click(button_clicked5)
        
        

        def rareplot(tipo_kits_1, tam_plot1, width_linea, linea_color, alfa_linea, size_point, alfa_point, markers_point, size_text_legend, size_title_legend, num_cols, family1, incluir_cuentas,
             label_X, ticklabels_X, label_Y, ticklabels_Y, family_axis, cambiar_sample, ajustar_ejes, Multiple_colorS, ver_num_muestra):
            
            def set_locus():
                mpl.rcParams.update(mpl.rcParamsDefault)
                AX1 = fig1.add_axes([0, 0, 1, 1])
                
                if tipo_kits_1 == 'Both kits':
                    NEW_RAREFACTION = RAREFACTION
                    sam_sacatter_counts = {i : RAREFACTION[i]['xdata'].max() for i in RAREFACTION}
                    val_max_x1 = max(set([j for i in RAREFACTION for j in RAREFACTION[i]['xdata']]))
                    val_max_y1 = max(set([j for i in RAREFACTION for j in RAREFACTION[i]['ydata']]))
                    
                    
                if tipo_kits_1 == 'DNPowerSoil':
                    NEW_RAREFACTION = {i : RAREFACTION[i] for i in RAREFACTION if i in DNPowerSoil}
                    sam_sacatter_counts = {i : RAREFACTION[i]['xdata'].max() for i in DNPowerSoil}
                    val_max_x1 = max(set([j for i in DNPowerSoil for j in RAREFACTION[i]['xdata']]))
                    val_max_y1 = max(set([j for i in DNPowerSoil for j in RAREFACTION[i]['ydata']]))
                    
                    
                if tipo_kits_1 == 'DNMicrobial':
                    NEW_RAREFACTION = {i : RAREFACTION[i] for i in RAREFACTION if i in DNMicrobial}
                    sam_sacatter_counts = {i : RAREFACTION[i]['xdata'].max() for i in DNMicrobial}
                    val_max_x1 = max(set([j for i in DNMicrobial for j in RAREFACTION[i]['xdata']]))
                    val_max_y1 = max(set([j for i in DNMicrobial for j in RAREFACTION[i]['ydata']]))
                    
                    
                AX1.clear()
                AX1 = fig1.add_axes([0, 0, 1, 1])
                
                
                AX1.tick_params(bottom=True, right=False, top=False, left=True, width = 2, length=4, color='gainsboro')
                
                if ajustar_ejes == 'True':
                    AX1.set_xlim(-val_max_x1*0.015, val_max_x1+(val_max_x1*0.03))
                    AX1.set_ylim(-val_max_y1*0.01, val_max_y1 + (val_max_y1 * 0.05))
                if ajustar_ejes == 'False':
                    AX1.set_xlim(-val_max_x*0.015, val_max_x+(val_max_x*0.03))
                    AX1.set_ylim(-val_max_y*0.01, val_max_y + (val_max_y * 0.05))
                
                
                fig1.set_size_inches((tam_plot1 * 1.5, tam_plot1))
                
                ### leyenda
                
                rampa_creada = list(dict.fromkeys([to_hex(k) for i in Multiple_colorS for k in plt.get_cmap(i)(np.arange(qualitative_colors[i])/qualitative_colors[i])]))
                
                cuenta_colors_tax = HBox([blanca, widgets.Label('MARKER COLORS='+str(len(rampa_creada))+', '+'SAMPLES='+str(len(KITS[tipo_kits_1])))])
                
                if len(KITS[tipo_kits_1]) <= len(rampa_creada):
                    rampa_creada = rampa_creada[0:len(KITS[tipo_kits_1])]
                else:
                    pass
                
                
                
                display(cuenta_colors_tax)
                
                cnorm = mpl.colors.Normalize(vmin=0, vmax=len(KITS[tipo_kits_1])-1)
                cpick = cm.ScalarMappable(norm=cnorm, cmap= ListedColormap(rampa_creada))
                ColormaP = [to_hex(cpick.to_rgba(e)) for e, i in enumerate(KITS[tipo_kits_1])]
                

                proxies = []
                labels = []
                for sam, col in zip(KITS[tipo_kits_1], ColormaP):
                    proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                             c=col, marker={'⏺':'o', '⏹':'s'}[markers_point], alpha = 1,
                                             markersize = size_point*12, url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[sam],
                                             markeredgecolor='none', linewidth = 0)
                    proxies.append(proxy)
                    if incluir_cuentas == 'True':
                        if cambiar_sample == 'Code1':
                            labels.append(sam+' ('+str(sam_sacatter_counts[sam])+')')
                        if cambiar_sample == 'Code2':
                            labels.append(name_code[sam]+' ('+str(sam_sacatter_counts[sam])+')')
                        if cambiar_sample == 'Code3':
                            labels.append(name_code2[sam]+' ('+str(sam_sacatter_counts[sam])+')')

                    if incluir_cuentas == 'False':
                        if cambiar_sample == 'Code1':
                            labels.append(sam)
                        if cambiar_sample == 'Code2':
                            labels.append(name_code[sam])
                        if cambiar_sample == 'Code3':
                            labels.append(name_code2[sam])
                
                
                constante = 15
                radio = size_point * constante # para normalizar la medida de 0 a 1
                factor1 = tam_plot1/(tam_plot1 * 1.5)
                factor2 = AX1.get_xlim()[1]/AX1.get_ylim()[1]
                factor3 = radio*factor1*factor2
                
                for co, j in zip(NEW_RAREFACTION, ColormaP): # ***********
                    xx = NEW_RAREFACTION[co]['xdata']
                    yy = NEW_RAREFACTION[co]['ydata']

                    line = mpl.lines.Line2D(xx.tolist(), yy.tolist(), linewidth=width_linea, alpha=alfa_linea, color = j, zorder=0)
                    AX1.add_line(line)

                    if markers_point == '⏺':
                        # elipse
                        AX1.add_patch(Ellipse((max(xx), max(yy)), factor3, radio, angle= 0, linewidth=0, facecolor=j, zorder=2, alpha = alfa_point,
                                              label = co, url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[co]))
                        
                    if markers_point == '⏹':
                        # rectangulo
                        AX1.add_patch(Rectangle((max(xx) - (factor3/2), max(yy) - (radio/2)), factor3, radio, angle= 0, linewidth=0, facecolor=j, alpha = alfa_point,
                                                zorder=2, label = co, url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[co]))
                        
                    if ver_num_muestra == 'True':
                        AX1.text(max(xx), max(yy), '   '+co.split('_')[-1], color = 'black', fontsize = 9, rotation = 0, ha='left', va="center", fontfamily = '3ds Light')
                        
                AX1.legend(proxies, labels, numpoints=1, title="$\\bf{Sample}$", title_fontsize = size_title_legend, loc=2, ncol = num_cols, facecolor = 'none',
                           fancybox=True, framealpha=1, shadow=False,
                           handlelength = 0.7, # ancho del marker
                           labelspacing = 0.6, # espacio entre filas
                           columnspacing = 1, # espacio entre colimnas
                           handletextpad=0.6, # espacio entre marker y label
                           borderpad = 0.5,
                           edgecolor="gainsboro", #frameon=False,
                           bbox_to_anchor=(1.01, 1.02),
                           prop={'family': family1, 'size': size_text_legend})
                
                AX1.set_xlabel('Contigs', fontsize=label_X, fontname=family_axis, weight="bold")
                
                AX1.set_ylabel('Species', fontsize=label_Y, fontname=family_axis, weight="bold")

                for labejex in AX1.xaxis.get_ticklabels():
                    labejex.set_fontsize(ticklabels_X)
                for labejey in AX1.yaxis.get_ticklabels():
                    labejey.set_fontsize(ticklabels_Y)

                for tickx in AX1.xaxis.get_major_ticks():
                    tickx.label.set_fontfamily([family_axis]) 
                for ticky in AX1.yaxis.get_major_ticks():
                    ticky.label.set_fontfamily([family_axis])

                AX1.xaxis.get_label().set_fontsize(label_X)
                AX1.yaxis.get_label().set_fontsize(label_Y)
                AX1.xaxis.get_label().set_fontfamily([family_axis])
                AX1.yaxis.get_label().set_fontfamily([family_axis])
                

                display(fig1)
                
            set_locus()
        
        
        out_rareplot = widgets.interactive_output(rareplot, {'tipo_kits_1':tipo_kits_1, 'tam_plot1':tam_plot1, 'width_linea':width_linea, 'linea_color':linea_color, 'alfa_linea':alfa_linea, 'size_point':size_point,
                                             'alfa_point':alfa_point, 'markers_point':markers_point,
                                             'size_title_legend':size_title_legend, 'size_text_legend':size_text_legend, 'num_cols':num_cols,
                                             'family1':family1, 'incluir_cuentas':incluir_cuentas,
                                             'label_X':label_X, 'ticklabels_X':ticklabels_X, 'label_Y':label_Y, 'ticklabels_Y':ticklabels_Y,
                                             'family_axis':family_axis, 'cambiar_sample':cambiar_sample, 'ajustar_ejes':ajustar_ejes, 'Multiple_colorS':Multiple_colorS, 'ver_num_muestra':ver_num_muestra})
        items_plot_1 = VBox([tam_plot1,
                             width_linea,
                             alfa_linea,
                             size_point,
                             alfa_point,
                             HBox([blanca, widgets.Label('Marker:'), markers_point, widgets.Label('Sample number:'), ver_num_muestra])])
        items_plot_1_box = Box(children=[items_plot_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_plot_1.children) * 30))+'px'))


        items_legend_1 = VBox([size_title_legend,
                               size_text_legend,
                               HBox([widgets.Label('Legend font:'), family1]),
                               HBox([blanca, widgets.Label('Columns:'), num_cols, blanca, HBox([widgets.Label('Contigs:'), incluir_cuentas])]),
                               HBox([blanca, widgets.Label('Rename:'), cambiar_sample])])
        items_legend_1_box = Box(children=[items_legend_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_legend_1.children) * 32))+'px'))


        items_axis_1 = VBox([label_X, ticklabels_X, label_Y, ticklabels_Y, HBox([blanca, widgets.Label('Axis font:'), family_axis])])
        items_axis_1_box = Box(children=[items_axis_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_axis_1.children) * 30))+'px'))
        
        
        items_save_1 = VBox([HBox([widgets.HTML('<font color = grey> <b style="font-size:0.7vw">SAVE CHART: </b>'), blanca, filename_plot]),
                             HBox([blanca, widgets.Label('Formats:'), png1, jpeg1, svg1, pdf1]),
                     HBox([blanca, widgets.Label('Chart parameters:'), filename, params])])
        items_save_1_box = Box(children=[items_save_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_save_1.children) * 34))+'px'))
        

        items_data = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'DATA']
        items_legend = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'BOX']
        items_axis = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'AXIS']
        items_save = [widgets.HTML('<font color = white> <b style="font-size:0.7vw">'+i+'</b>') for i in 'S']
        
        
        RAREFY = HBox([VBox([widgets.HTML('<font color = #1976d2> <b style="font-size:0.8vw">KITS</b>'), tipo_kit_1_box,
                             VBox([widgets.HTML('<font color = gray> <h style="font-size:0.55vw">Adjust axes</h>'), ajustar_ejes])]), blanca,
                       VBox([#widgets.HTML('<i class="fa fa-cog fa-2x fa-fw"></i>  <b style="font-size:1vw">PLOT SETTINGS</b>'),
                                             #widgets.HTML('<i class="fa fa-cog fa-spin fa-2x fa-fw"></i> <b style="font-size:1vw">PLOT SETTINGS</b>'),
                                             widgets.HTML('<font color = grey> <i class="fa fa-cog fa-2x fa-fw"></i> <b style="font-size:0.8vw">PLOT SETTINGS</b>'),
                                             #widgets.HTML('<b style="font-size:0.8vw">PLOT SETTINGS</b>'),
                                             
                                             #widgets.Label('DATA'),
                                             HBox([Box(children =[VBox(items_data)]), items_plot_1_box]),
                           
                                             #widgets.Label('AXIS'),
                                             HBox([Box(children =[VBox(items_axis)]), items_axis_1_box]),
                           
                                             #widgets.Label('LEGEND'),
                                             HBox([Box(children =[VBox(items_legend)]), items_legend_1_box]),
                                             
                                             #widgets.Label('SAVE'),
                                             HBox([Box(children =[VBox(items_save)]), items_save_1_box])
                       ]),
                       VBox([widgets.HTML('<font color = grey> <b style="font-size:0.6vw">MARKER COLORS</b>'), Multiple_colorS, HBox([blanca, blanca, ayuda3])]), blanca, VBox([out_rareplot])])

        display(RAREFY, output1)
        

RARE_button.on_click(button_clicked)





RARE_ANALYSIS = VBox([HBox([RARE_button, estatico_box]), RARE_output])








# # Alpha indices




help1_button = widgets.Button(description="Help", icon = 'fa-question-circle', layout=Layout(width='65px'))
help1_button.style.button_color = 'pink'
help1_button.style.font_weight = 'bold'
help1_button_output = widgets.Output()

import tkinter as tk
from tkinter import *
    
def button_clicked(b):
    with help1_button_output:
        clear_output(True)
        
        #---------------------------
        #import tkinter as tk
        #from tkinter import messagebox
        #import ctypes
        #ctypes.windll.shcore.SetProcessDpiAwareness(1)

        #root = tk.Tk()
        #root.overrideredirect(1)
        #root.withdraw()
        #root.attributes('-topmost', 1)
        #messagebox.showinfo("Information about phylogenetic diversity analysis", "For this analysis, the most representative ASVs of each species were taken (with more readings), "\
        #                    "from these a fasta file is obtained that is aligned using Clustalo.\n\nClustalo v1.2.2\n\n"\
        #                    "Fast, scalable generation of high‐quality protein multiple sequence alignments using Clustal Omega\n\ndoi: 10.1038/msb.2011.75")
        #root.destroy()
        #-----------------------------------------
        
        import webbrowser
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)

        newwin0 = Tk()
        newwin0.title("Information about analysis")

        newwin0.geometry("515x480")
        newwin0.configure(background='white')
        newwin0.resizable(0, 0)
        newwin0.attributes('-topmost', 1)

        c00 = Label(newwin0, text="    ", bg = 'white')
        c00.grid(column=0, row=0)

        text2 = Text(newwin0, height=40, width=80, bg = 'white', bd = 0, font=('Arial', 7))
        scroll = Scrollbar(newwin0, command=text2.yview)
        text2.configure(yscrollcommand=scroll.set)

        text2.tag_configure('big', font=('Arial', 8, 'bold'))
        # titulos

        # contenidos
        text2.tag_configure('color', foreground='red', font=('Arial', 7, 'bold'))
        text2.tag_configure('color2', foreground='darkorange', font=('Arial', 7, 'bold'))

        text2.tag_configure('text', font=('Arial', 7))
        text2.tag_configure('text2', font=('Arial', 7, 'bold'))

        hyperlink = HyperlinkManager(text2)

        #------------------------------------------------------------
        text2.insert(END,'\nPhylogenetic diversity\n\n', 'big')
        
        text2.insert(END, "This analysis requires the installation of scikit-bio v0.5.6\n", 'text')
        def click1():
            webbrowser.open_new(r"http://scikit-bio.org/")
        text2.insert(END, "http://scikit-bio.org/\n\n", hyperlink.add(click1))
        
        
        text2.insert(END, "For this analysis the ASVs were collapsed to species taking those with\n"                     "the most readings, from these a fasta file was built\n"                     "which was aligned using Clustalo.\n\n", 'text')
        

        text2.insert(END, "\nThis approach was adapted from scikit-bio and IAB:\n")
        def click1():
            webbrowser.open_new(r"http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.faith_pd.html")
        text2.insert(END, "http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.faith_pd.html\n", hyperlink.add(click1))
        def click1():
            webbrowser.open_new(r"http://readiab.org/book/0.1.3/3/1")
        text2.insert(END, "http://readiab.org/book/0.1.3/3/1\n", hyperlink.add(click1))


        text2.insert(END, "\n\nDownload Clustalo v1.2.2\n", 'text')
        def click1():
            webbrowser.open_new(r"http://www.clustal.org/omega/")
        text2.insert(END, "http://www.clustal.org/omega/", hyperlink.add(click1))

        text2.insert(END, "\n\nFor more information about Clustalo visit:\n", 'text')

        text2.insert(END, "\nFast, scalable generation of high quality protein multiple sequence alignments \n"                     "using Clustal Omega\n", 'text')

        def click1():
            webbrowser.open_new(r"https://doi.org/10.1038/msb.2011.75")
        text2.insert(END, "https://doi.org/10.1038/msb.2011.75", hyperlink.add(click1))
        
        text2.insert(END, "\n\n-------------------------------------------------------------------------------------------------", 'text')
        
        text2.insert(END, "\n\nYou can see the file called ASVs_Phylogenetic_tree_VIEW_Species.txt\n"                     "deposited in the clustering folder.\n", 'text')
        
        def click1():
            webbrowser.open_new(r"https://itol.embl.de/upload.cgi")
        text2.insert(END, "https://itol.embl.de/upload.cgi", hyperlink.add(click1))


        text2.grid(column=1, row=0)
        scroll.grid(column=1, row=0, sticky = 'ESN')
        text2.configure(state=DISABLED)
        newwin0.mainloop()
        
         
help1_button.on_click(button_clicked)
ayuda1 = HBox([help1_button, help1_button_output])





alfa_filo_button = widgets.Button(description="RUN: Philogenetic diversity", icon = 'fa-play', layout=Layout(width='230px'))
alfa_filo_button.style.button_color = 'beige'
alfa_filo_button.style.font_weight = 'bold'
alfa_filo_output = widgets.Output()

def button_clicked(b):
    
    from io import StringIO
    from skbio.tree import TreeNode
        
    NCBI_RDP_SILVA_SUMMARY = pd.read_csv('tablas/ASVs_NCBI_RDP_SILVA_SUMMARY.txt', sep = '\t')
    
    rarefaction_collapse_specie = pd.pivot_table(NCBI_RDP_SILVA_SUMMARY[['Species'] + name_sample], values = name_sample, index = ['Species'], aggfunc = sum).reset_index()
    
    with alfa_filo_output:
        clear_output(True)
        
        progreso = widgets.HTML()
        progress = HBox([widgets.HTML('<font color = black> <i class="fa fa-spinner fa-pulse fa-2x fa-fw"></i> </font>'), progreso])
        display(Box(children=[progress], layout=Layout(border='1px solid limegreen', width='585px', height='32px')))

        
        progreso.value = '<font color = limegreen> <b style="font-size:0.5vw">Process: </b> <font color = black> <b style="font-size:0.5vw"> Processing Data </b>'
        
        """
        http://readiab.org/book/0.1.3/
        http://readiab.org/book/0.1.3/3/1#3

        For comparing phylogenetic diversity among communities, go ahead and use synthesis phylogenies.
        PD calculates the sum of the branch lengths of all species present in an assemblage.
        We didnotinclude the root of the phylogeny when calculating PD.

        Phylogenetic diversity (Faith’s PD) uses phylogenetic distance to calculate the diversity of a given sample.
        """
        def get_observed_nodes(tree, table, sample_id, columna, verbose=False):

            observed_otus = table[[sample_id, columna]][table[sample_id] > 0][columna].tolist()

            #observed_otus = [obs_id for obs_id in table.index if table[sample_id][obs_id] > 0]

            observed_nodes = set()
            # iterate over the observed OTUs
            for otu in observed_otus:
                t = tree.find(otu)
                observed_nodes.add(t)
                if verbose:
                    print(t.name, t.length, end=' ')
                for internal_node in t.ancestors():
                    if internal_node.length is None:
                        # we've hit the root
                        if verbose:
                            pass
                    else:
                        if verbose and internal_node not in observed_nodes:
                            print(internal_node.length, end=' ')
                        observed_nodes.add(internal_node)
            return observed_nodes

        def phylogenetic_diversity(tree, table, sample_id, columna, verbose=False):
            observed_nodes = get_observed_nodes(tree, table, sample_id, columna, verbose=verbose)
            result = sum(o.length for o in observed_nodes)
            return result

        fasta = open('clustering/ASVs.fasta')
        fastas = fasta.read()
        fasta.close()
        def find_seq(iden = ''):
            for i in fastas.split('>')[1:]:
                i = i.rstrip()
                identificador = re.search('\w+', i).group()
                if iden == identificador:
                    secuencia = re.sub(identificador, '', i)
                    secuencia = re.sub('\n', '', secuencia)
                    secuencia = re.sub('[*]', '', secuencia)
                    return secuencia
                    break

        asv_species_dict = {}
        asv_species_lista = []
        for i in NCBI_RDP_SILVA_SUMMARY.Species.drop_duplicates():
            w = NCBI_RDP_SILVA_SUMMARY[NCBI_RDP_SILVA_SUMMARY.Species == i]
            sumas = np.sum(w.iloc[:, 3:].values, axis = 1)
            w['sum'] = sumas
            w = w.sort_values(by =['sum'],ascending=False).reset_index(drop=True)
            ef = w['Entry'].tolist()[0]
            asv_species_dict[ef] = i
            asv_species_lista.append([ef, i])
        unicos = DataFrame(asv_species_lista, columns = ['Entry', 'Species'])

        MERGE_UNICOS_SPECIES = unicos.merge(rarefaction_collapse_specie, on = 'Species', how = 'left')
        
        progreso.value = '<font color = limegreen> <b style="font-size:0.5vw">Process: </b> <font color = black> <b style="font-size:0.5vw"> ASVs_Phylogenetic_Diversity.fasta </b>'
        
        with open('clustering/ASVs_Phylogenetic_Diversity.fasta', 'w') as fq:
            for e, i in enumerate(asv_species_dict):
                if e == (len(asv_species_dict)-1):
                    fq.write('>'+i+'\n'+find_seq(iden = i))
                else:
                    fq.write('>'+i+'\n'+find_seq(iden = i)+'\n')
                    
        with open('tablas/asv_species_dict.json', 'w') as fp:
                json.dump(asv_species_dict, fp)
                
                
        
        progreso.value = '<font color = limegreen> <b style="font-size:0.5vw">Process: </b> <font color = black> <b style="font-size:0.5vw"> Building Alignments  </b>'

        """
        Alineamiento multiple de secuencias usando clustalo, forzado a sobreescribir archivos
        """
        #----------------
        info10 = subprocess.Popen('binarios/clustalo.exe --infile clustering/ASVs_Phylogenetic_Diversity.fasta --threads 8 --MAC-RAM 8000 --verbose --guidetree-out clustering/ASVs_Phylogenetic_Diversity_tree.txt --outfmt clustal --force --resno --outfile clustering/ASVs_Phylogenetic_Diversity_alig.txt --output-order tree-order --seqtype dna',
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        info100 = info10.communicate()[0].decode()
        

        arbol1 = open('clustering/ASVs_Phylogenetic_Diversity_tree.txt', 'r')
        arbol = arbol1.read()
        arbol1.close()
        arbol = ''.join(arbol.split('\n'))
        arbol = re.sub(';$', 'root;', arbol)

        # este arbol filogenetico lo obtuve desde clustal omega
        newick_tree = StringIO(arbol)

        tree = TreeNode.read(newick_tree)
        tree = tree.root_at_midpoint()
        
        #with open('arbol.txt', 'w') as fq:
        #    fq.write(tree.ascii_art())

        PD = []
        for i in name_sample:
            pd_A = phylogenetic_diversity(tree, MERGE_UNICOS_SPECIES, i, 'Entry', verbose=False)
            PD.append([pd_A])
        PD = [['Phylogenetic Diversity']] + PD
        PD = DataFrame(np.array(PD).T, columns = ['Index'] + name_sample)
        PD.to_csv('tablas/ASVs_PD.txt',index = None, sep = '\t')

        clear_output(True)
        progreso.value = '<font color = limegreen> <b style="font-size:0.5vw">Process: </b> <font color = black> <b style="font-size:0.5vw"> Finished </b>'
        progress = HBox([widgets.HTML('<font color = black> <i class="fa fa-spinner fa-2x fa-fw"></i> </font>'), progreso])
        display(Box(children=[progress], layout=Layout(border='1px solid limegreen', width='585px', height='32px')))

alfa_filo_button.on_click(button_clicked)





view_tree_button = widgets.Button(description="View phylogenetic tree in plain text ", icon = 'fa-eye', layout=Layout(width='280px'))
view_tree_button.style.button_color = 'aquamarine'
view_tree_button.style.font_weight = 'bold'
view_tree_output = widgets.Output()

def button_clicked(b):
    from io import StringIO
    from skbio.tree import TreeNode
    
    with open('tablas/asv_species_dict.json', 'r') as fp:
        asv_species_dict = json.load(fp)
    
    with view_tree_output:
        clear_output(True)

        arbol11 = open('clustering/ASVs_Phylogenetic_Diversity_tree.txt', 'r')
        arb = ''
        for line in arbol11:
            line = line.rstrip()
            if re.search('ASV\d+', line):
                arb += re.sub('ASV\d+', asv_species_dict[re.search('ASV\d+', line).group()], line)
            else:
                arb += line
        arb = re.sub(';$', 'root;', arb)
        
        with open('clustering/ASVs_Phylogenetic_tree_VIEW_Species.txt', 'w') as fq:
            fq.write(arb)
            
        newick_tree2 = StringIO(arb)

        tree2 = TreeNode.read(newick_tree2)
        tree2 = tree2.root_at_midpoint()

        with open('clustering/ASVs_Phylogenetic_tree_VIEW.txt', 'w') as fq:
            fq.write(tree2.ascii_art())

        os.system('start notepad.exe clustering/ASVs_Phylogenetic_tree_VIEW.txt')

view_tree_button.on_click(button_clicked)





diver_filogen = VBox([HBox([alfa_filo_button, ayuda1, view_tree_button, view_tree_output]), alfa_filo_output])





help2_button = widgets.Button(description="Help", icon = 'fa-question-circle', layout=Layout(width='65px'))
help2_button.style.button_color = 'pink'
help2_button.style.font_weight = 'bold'
help2_button_output = widgets.Output()

import tkinter as tk
from tkinter import *
    
def button_clicked(b):
    with help2_button_output:
        clear_output(True)
        
        import webbrowser
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)

        newwin0 = Tk()
        newwin0.title("Information about analysis")

        newwin0.geometry("515x280")
        newwin0.configure(background='white')
        newwin0.resizable(0, 0)
        newwin0.attributes('-topmost', 1)

        c00 = Label(newwin0, text="    ", bg = 'white')
        c00.grid(column=0, row=0)

        text2 = Text(newwin0, height=40, width=80, bg = 'white', bd = 0, font=('Arial', 7))
        scroll = Scrollbar(newwin0, command=text2.yview)
        text2.configure(yscrollcommand=scroll.set)

        text2.tag_configure('big', font=('Arial', 8, 'bold'))
        # titulos

        # contenidos
        text2.tag_configure('color', foreground='red', font=('Arial', 7, 'bold'))
        text2.tag_configure('color2', foreground='darkorange', font=('Arial', 7, 'bold'))

        text2.tag_configure('text', font=('Arial', 7))
        text2.tag_configure('text2', font=('Arial',7, 'bold'))

        hyperlink = HyperlinkManager(text2)

        #------------------------------------------------------------
        text2.insert(END,'\nDiversity calculations\n\n', 'big')
        
        text2.insert(END, '\nThis option provides functionality for analyzing biological diversity.\n\n', 'text')
        
        text2.insert(END, "\nThis analysis was adapted from the metrics used in the PAST data analysis\n"                     " program v4.03\n\n", 'text')
        
        def click1():
            webbrowser.open_new(r"http://priede.bf.lu.lv/ftp/pub/TIS/datu_analiize/PAST/2.17c/doc1.html")
        text2.insert(END, "http://priede.bf.lu.lv/ftp/pub/TIS/datu_analiize/PAST/2.17c/doc1.html\n\n", hyperlink.add(click1))
        
        text2.insert(END, '\nDownload PAST v4.XX:\n', 'text')
        
        def click1():
            webbrowser.open_new(r"https://www.nhm.uio.no/english/research/infrastructure/past/")
        text2.insert(END, "https://www.nhm.uio.no/english/research/infrastructure/past/\n\n", hyperlink.add(click1))
        
        text2.grid(column=1, row=0)
        scroll.grid(column=1, row=0, sticky = 'ESN')
        text2.configure(state=DISABLED)
        newwin0.mainloop()
        
         
help2_button.on_click(button_clicked)
ayuda2 = HBox([help2_button, help2_button_output])





index_names = ['Taxa_S','Individuals','Dominance_D','Simpson_1-D','Shannon_H','Evenness_e^H/S',
               'Brillouin','Menhinick','Margalef','Equitability_J','Berger-Parker','Chao-1']
"""
Diversity index as well as in PAST software 4.03
"""
def TAXA_S(entrada = int()):
    e = entrada
    return e
def INDIVIDUALS(entrada = list()):
    e = sum(entrada)
    return e
def DOMINANCE_D(entrada = list(), individuos = int()):
    din = []
    for i in entrada:
        din.append((i/individuos)**2)
    return sum(din)
def SIMPSON_1_D(entrada = float()):
    e = 1 - entrada
    return e
def SHANON_H(entrada = list(), individuos = int()):
    din = []
    for i in entrada:
        din.append((i/individuos) * np.log(i/individuos))
    return round(-sum(din), 3)
def EVENNESS_E_HS(shanon = float(), taxa = int()):
    e = (math.e**shanon)/taxa
    return e
def BRILLOUIN(entrada = list(), individuos = int()):
    din = []
    for i in entrada:
        din.append(math.log(math.factorial(i)))
    e = (math.log(math.factorial(individuos)) - sum(din)) / individuos
    return e
def MENHINICK(taxa = int(), individuos = int()):
    e = taxa / np.sqrt(individuos)
    return e
def MARGALEF(taxa = int(), individuos = int()):
    e = (taxa - 1) / np.log(individuos)
    return e
def EQUITABILITY_J(shanon = float(), taxa = int()):
    e = shanon / np.log(taxa)
    return e
#def FISHER_ALPHA(taxa = int(), individuos = int()):
#    from scipy.optimize import minimize_scalar
#    
#    def f(alpha):
#        return (alpha * np.log(1 + (individuos / alpha)) - taxa) ** 2
#    
#    e = minimize_scalar(f).x
#    return e
def BERGER_PARKER(entrada = list(), individuos = int()):
    e = round(max(entrada) / individuos, 4)
    return e
def CHAO_1(entrada = list(), taxa = int()):
    F1 = len([i for i in entrada if i == 1])
    F2 = len([i for i in entrada if i == 2])
    e = taxa + F1 * (F1 - 1) / (2 * (F2 + 1))
    return int(e)











indices_button = widgets.Button(description="RUN: Diversity calculations", icon = 'fa-play', layout=Layout(width='230px'))
indices_button.style.button_color = 'beige'
indices_button.style.font_weight = 'bold'
indices_button_output = widgets.Output()


def button_clicked(b):
    
    with open('tablas/ASVs_data_for_rarefaction.json', 'r') as fp:
        data_for_rarefaction = json.load(fp)
        
    import time
    
    with indices_button_output:
        clear_output(True)
        
        progreso = widgets.HTML()
        progress = HBox([widgets.HTML('<font color = black> <i class="fa fa-spinner fa-pulse fa-2x fa-fw"></i> </font>'), progreso])
        display(Box(children=[progress], layout=Layout(border='1px solid limegreen', width='585px', height='32px')))

        indices_diversidad = {}
        for sam in data_for_rarefaction:
            taxa_s = TAXA_S(entrada = data_for_rarefaction[sam][0]['No_species'])
            individuals = INDIVIDUALS(entrada = data_for_rarefaction[sam][1]['Counts'])
            dominance_d = DOMINANCE_D(entrada = data_for_rarefaction[sam][1]['Counts'], individuos = individuals)
            simpson_1_d = SIMPSON_1_D(entrada = dominance_d)
            shanon_h = SHANON_H(data_for_rarefaction[sam][1]['Counts'], individuos = individuals)
            evenness_e_hs = EVENNESS_E_HS(shanon = shanon_h, taxa = taxa_s)
            brillouin = BRILLOUIN(data_for_rarefaction[sam][1]['Counts'], individuos = individuals)
            menhinick = MENHINICK(taxa = taxa_s, individuos = individuals)
            margalef = MARGALEF(taxa = taxa_s, individuos = individuals)
            equitability_j = EQUITABILITY_J(shanon = shanon_h, taxa = taxa_s)
            #fisher_alpha = FISHER_ALPHA(taxa = taxa_s, individuos = individuals)
            berger_parker = BERGER_PARKER(entrada = data_for_rarefaction[sam][1]['Counts'], individuos = individuals)
            chao_1 = CHAO_1(entrada = data_for_rarefaction[sam][1]['Counts'], taxa = taxa_s)
            indices_diversidad[sam] = [taxa_s,individuals,dominance_d,simpson_1_d,shanon_h,evenness_e_hs,brillouin,
                                      menhinick,margalef,equitability_j,berger_parker,chao_1]
            
            
            progreso.value = '<font color = limegreen> <b style="font-size:0.5vw">Processed samples : </b> <font color = black> <b style="font-size:0.5vw"> '+sam+' </b>'
            time.sleep(0.2)
        
        clear_output(True)
        progreso.value = '<font color = limegreen> <b style="font-size:0.5vw">Processed samples : </b> <font color = black> <b style="font-size:0.5vw"> '+sam+' </b>'
        progress = HBox([widgets.HTML('<font color = black> <i class="fa fa-spinner fa-2x fa-fw"></i> </font>'), progreso])
        display(Box(children=[progress], layout=Layout(border='1px solid limegreen', width='585px', height='32px')))

        Indices = pd.DataFrame(data=indices_diversidad)
        Indices.insert(loc = 0, column='Index', value=index_names)
        Indices.to_csv('tablas/ASVs_Indices.txt',index = None, sep = '\t')

        
indices_button.on_click(button_clicked)





metricas_diversidad = VBox([HBox([indices_button, ayuda2]), indices_button_output])





PHY_DIV = HBox([diver_filogen, blanca, metricas_diversidad])

















INDICES_button = widgets.Button(description="PROCESS AND VISUALIZE", icon = 'fa-eye', layout=Layout(width='590px'))
INDICES_button.style.button_color = 'gainsboro' #'deepskyblue'
INDICES_button.style.font_weight = 'bold'
INDICES_output = widgets.Output()

def button_clicked(b):
    
    import os.path
    from matplotlib.colors import ListedColormap # funcion para crear un objeto <matplotlib.colors.ListedColormap> a partir de una lista de colores personalizados
    
    if os.path.isfile('tablas/ASVs_PD.txt'):
        PD  = pd.read_csv('tablas/ASVs_PD.txt', sep = '\t')
        Indices = pd.read_csv('tablas/ASVs_Indices.txt', sep = '\t')
        Indices = pd.concat([Indices, PD]).reset_index(drop = True)
    else:
        Indices = pd.read_csv('tablas/ASVs_Indices.txt', sep = '\t')
    
    with INDICES_output:
        clear_output(True)
        
        mpl.rcParams.update(mpl.rcParamsDefault)
        fig2 = plt.figure(figsize=(6, 4))

        AX2 = fig2.add_axes([0, 0, 1, 1])

        AX2.set_facecolor('none')
        
        AX2.tick_params(bottom=True, right=False, top=False, left=True, width = 2, length=4, color='gainsboro')
        plt.gca().tick_params(which='major', width = 2, length=4, color='gainsboro')
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_color('gainsboro')
        plt.gca().spines['bottom'].set_color('gainsboro')
        plt.gca().spines['right'].set_color(None)
        plt.gca().spines['top'].set_color(None)
        
        
        plt.close()
        
        tipo_kits_1 = widgets.ToggleButtons(options= ['Both kits'] + metadata[VARIABLE_KIT].unique().tolist(), value = 'Both kits', button_style = 'primary')
        tipo_kits_1.style.button_width = '170px'
        tipo_kits_1.style.font_weight = 'bold'
        tipo_kit_1_box = Box(children=[VBox([tipo_kits_1])], layout=Layout(border='1px solid #1976d2', width='180px', height='95px'))
        
        tipo_indice = widgets.ToggleButtons(options= list(reversed(list(dict(OrderedDict(Counter({i : len(i) for i in Indices.Index.tolist()}).most_common())).keys()))), value = 'Chao-1', button_style = '')
        tipo_indice.style.button_width = '170px'
        tipo_indice_box = Box(children=[VBox([tipo_indice])], layout=Layout(border='1px solid #1976d2', width='180px', height='380px'))
        
        Multiple_colorS = widgets.SelectMultiple(options=sorted(list(qualitative_colors.keys())), value=sorted(list(qualitative_colors.keys()))[0:5], disabled=False,
                       layout=Layout(width='110px', height='250px'))
        
        fig_size = []
        n = 2
        for i in range(16):
            fig_size.append(round(n, 2))
            n += 0.2
    
        tam_plot1 = widgets.SelectionSlider(options=fig_size,value=3,disabled=False,
                                              description = 'Chart size:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        ## lineas
        width_linea = widgets.SelectionSlider(options=width_line,value=3,disabled=False,
                                              description = 'Line width:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                             layout=Layout(width='400px', height='25px'))
        linea_color = widgets.ColorPicker(concise=True, value='black', disabled=False, layout = Layout(width='35px', height='25px'))
        
        
        ## lineas
        width_linea = widgets.SelectionSlider(options=width_line,value=2,disabled=False,
                                              description = 'Line width:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                             layout=Layout(width='400px', height='25px'))
        
        linea_color = widgets.ColorPicker(concise=True, value='#c0c0c0', disabled=False, layout = Layout(width='35px', height='25px'))
        
        text_size = widgets.SelectionSlider(options=range(0, 31),value=10,disabled=False,
                                              description = 'Text size:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        
        text_color = widgets.ColorPicker(concise=True, value='#000000', disabled=False, layout = Layout(width='35px', height='25px'))
        
        
        text_alfa = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=1,disabled=False,
                                              description = 'Text alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        alfa_linea = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=1,disabled=False,
                                              description = 'Line alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        # puntos

        size_point = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.5,disabled=False,
                                              description = 'Marker size:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        #point_color = widgets.ColorPicker(concise=True, value='lime', disabled=False, layout = Layout(width='35px', height='25px'))
        
        #formas = {'Circle': 'o', 'Square': 'o'}
        #markers_point = widgets.ToggleButtons(options=list(formas.keys()), value = 'Square')
        
        #markers_point = widgets.ToggleButtons(options=['⚫', '⬛'], value = '⚫')
        #markers_point.style.button_width = '40px'
        
        markers_point = widgets.ToggleButtons(options=['⏺', '⏹'], value = '⏺')
        markers_point.style.button_width = '31px'
        
        alfa_marker = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=1,disabled=False,
                                              description = 'Marker alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        
        family = sorted(['Liberation Serif','Microsoft Sans Serif','Open Sans','Times New Roman','3ds Light','Calibri','Comic Sans MS',
                          'Arial','Courier New','Microsoft Yi Baiti','Lucida Console'])
        family1 = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                           layout = Layout(width='290px', height='25px'))
        
        ### ejes
        label_X = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'X label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_X = widgets.SelectionSlider(options=range(5, 31),value=10,disabled=False,
                                              description = 'X tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))

        label_Y = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'Y label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_Y = widgets.SelectionSlider(options=range(5, 31),value=10,disabled=False,
                                              description = 'Y tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))

        family_axis = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                   layout = Layout(width='290px', height='25px'))
        
        cambiar_sample = widgets.ToggleButtons(options=['Code1', 'Code2', 'Code3'], value = 'Code1')
        cambiar_sample.style.button_width = '70px'
        
        ajustar_ejes = widgets.ToggleButtons(options=['True', 'False'], value = 'False')
        ajustar_ejes.style.button_width = '55px'
        
        #--------------------------
        
        filename_plot = widgets.Text(value='', placeholder='Chart name', description='', disabled=False, layout = Layout(width='270px', height='25px'))
        
        ancho_plot_save = str(71)
         
        png1 = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
        png1.style.button_color = 'gold'
        output1 = widgets.Output()
        def button_clicked1(b):
            with output1:
                clear_output(True)
                if filename_plot.value is '':
                    nombre_grafico = 'Indices_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                if filename_plot.value is not '':
                    nombre_grafico = filename_plot.value
                fig2.savefig('plots_asv/indices/'+nombre_grafico+'.png', dpi = 900, bbox_inches= 'tight')
                
        png1.on_click(button_clicked1)
        #----
        jpeg1 = widgets.Button(description="JPEG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
        jpeg1.style.button_color = 'gold'
        output2 = widgets.Output()
        def button_clicked2(b):
            with output2:
                clear_output(True)
                if filename_plot.value is '':
                    nombre_grafico = 'Indices_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                if filename_plot.value is not '':
                    nombre_grafico = filename_plot.value
                fig2.savefig('plots_asv/indices/'+nombre_grafico+'.jpeg', dpi = 900, bbox_inches= 'tight')
                
        jpeg1.on_click(button_clicked2)
        #----
        svg1 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
        svg1.style.button_color = 'gold'
        output3 = widgets.Output()
        def button_clicked3(b):
            with output3:
                clear_output(True)
                if filename_plot.value is '':
                    nombre_grafico = 'Indices_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                if filename_plot.value is not '':
                    nombre_grafico = filename_plot.value
                fig2.savefig('plots_asv/indices/'+nombre_grafico+'.svg', dpi = 900, bbox_inches= 'tight')
                
        svg1.on_click(button_clicked3)
        #----
        pdf1 = widgets.Button(description="PDF", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
        pdf1.style.button_color = 'gold'
        output4 = widgets.Output()
        def button_clicked4(b):
            with output4:
                clear_output(True)
                if filename_plot.value is '':
                    nombre_grafico = 'Indices_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                if filename_plot.value is not '':
                    nombre_grafico = filename_plot.value
                fig2.savefig('plots_asv/indices/'+nombre_grafico+'.pdf', dpi = 900, bbox_inches= 'tight')
                
        pdf1.on_click(button_clicked4)
        #----
        filename = widgets.Text(value='', placeholder='File name.txt', description='', disabled=False, layout = Layout(width='195px', height='25px'))
        
        params = widgets.Button(description="SAVE", icon = 'fa-save', layout=Layout(width='70px'))
        params.style.button_color = 'cyan'
        output5 = widgets.Output()
        def button_clicked5(b):
            with output5:
                clear_output(True)
                
                if filename.value is '':
                    nombre_archivo = 'Indices_chart_params_ASVs.txt'
                if filename.value is not '':
                    nombre_archivo = filename.value
                
                with open('plots_asv/indices/'+nombre_archivo, 'w') as fq:
                    fq.write('#Saved parameters\n')
                    fq.write('Clustering:ASVs\n' )
                    """
                    diccionario para guardar en un archivo los parametros usados para volver a reproducirlos posteriormente
                    """
                    parametros = {'Kit:':tipo_kits_1.value,
                                  'Index:':tipo_indice.value,
                                  'Chart size:':tam_plot1.value,
                                  'Line width:':width_linea.value,
                                  'Line alpha:':alfa_linea.value,
                                  'Text size:':text_size.value,
                                  'Text font:':family1.value,
                                  'Text alpha:':text_alfa.value,
                                  'Marker colors:':', '.join(list(Multiple_colorS.value)),
                                  'Marker size:':size_point.value,
                                  'Marker:':{'⏺':'circle', '⏹':'square'}[markers_point.value],
                                  'Line color:':linea_color.value+' | '+'('+', '.join([r+':'+str(int(255*i)) for i, r in zip(to_rgb(linea_color.value), ['RED', 'GREEN', 'BLUE'])])+')',
                                  'Text color:':text_color.value+' | '+'('+', '.join([r+':'+str(int(255*i)) for i, r in zip(to_rgb(text_color.value), ['RED', 'GREEN', 'BLUE'])])+')',
                                  'X label:':label_X.value,
                                  'X tick labels:':ticklabels_X.value,
                                  'Y label:':label_Y.value,
                                  'Y tick labels:':ticklabels_Y.value,
                                  'Axis font:':family_axis.value,
                                  'Rename sample:':cambiar_sample.value}
                    for w in parametros:
                        fq.write(w+str(parametros[w])+'\n')
        
        params.on_click(button_clicked5)
        #------------------------------------
        
        def updata_plot(tipo_kits_1, tipo_indice, tam_plot1, width_linea, linea_color, text_size, text_color, text_alfa, alfa_linea, size_point, markers_point, alfa_marker, family1,
                        label_X, ticklabels_X, label_Y, ticklabels_Y, family_axis, cambiar_sample, Multiple_colorS, ajustar_ejes):
            def set_locus():
                
                mpl.rcParams.update(mpl.rcParamsDefault)
                AX2 = fig2.add_axes([0, 0, 1, 1])
                
                if tipo_kits_1 == 'Both kits':
                    fig_per = 2.2
                    ejeY = Indices[Indices.Index == tipo_indice][KITS['Both kits']].values[0]
                    ejeX = Indices[Indices.Index == tipo_indice][KITS['Both kits']].columns
                    val_max_y = max(ejeY)
                    val_min_y = min(ejeY)
                    aumentoX = 0.02
                if tipo_kits_1 == 'DNPowerSoil':
                    fig_per = 1.1
                    ejeY = Indices[Indices.Index == tipo_indice][KITS['DNPowerSoil']].values[0]
                    ejeX = Indices[Indices.Index == tipo_indice][KITS['DNPowerSoil']].columns
                    val_max_y = max(ejeY)
                    val_min_y = min(ejeY)
                    aumentoX = 0.055
                if tipo_kits_1 == 'DNMicrobial':
                    fig_per = 1.1
                    ejeY = Indices[Indices.Index == tipo_indice][KITS['DNMicrobial']].values[0]
                    ejeX = Indices[Indices.Index == tipo_indice][KITS['DNMicrobial']].columns
                    val_max_y = max(ejeY)
                    val_min_y = min(ejeY)
                    aumentoX = 0.055
                
                
                
                rampa_creada = list(dict.fromkeys([to_hex(k) for i in Multiple_colorS for k in plt.get_cmap(i)(np.arange(qualitative_colors[i])/qualitative_colors[i])]))
                
                cuenta_colors_tax = HBox([blanca, widgets.Label('MARKER COLORS='+str(len(rampa_creada))+', '+'SAMPLES='+str(len(KITS[tipo_kits_1])))])
                
                if len(KITS[tipo_kits_1]) <= len(rampa_creada):
                    rampa_creada = rampa_creada[0:len(KITS[tipo_kits_1])]
                else:
                    pass
                
                display(cuenta_colors_tax)
                
                cnorm = mpl.colors.Normalize(vmin=0, vmax=len(KITS[tipo_kits_1])-1)
                cpick = cm.ScalarMappable(norm=cnorm, cmap= ListedColormap(rampa_creada))
                ColormaP = [to_hex(cpick.to_rgba(e)) for e, i in enumerate(KITS[tipo_kits_1])]
                
                
                
                fig2.set_size_inches((tam_plot1 * fig_per, tam_plot1)) # size plot
                
                AX2.clear()
                AX2 = fig2.add_axes([0, 0, 1, 1])
                
                AX2.set_xlim(-len(ejeY)*aumentoX, (len(ejeY)-1)+((len(ejeY)-1)*aumentoX))
                
                if ajustar_ejes == 'True':
                    AX2.set_ylim(val_min_y - (val_max_y * 0.05), val_max_y + (val_max_y * 0.05)) # AX2.set_ylim(-val_max_y*0.015, val_max_y + (val_max_y * 0.05))
                if ajustar_ejes == 'False':
                    AX2.set_ylim(0, val_max_y + (val_max_y * 0.05))
                
                
                line = mpl.lines.Line2D([e for e, i in enumerate(ejeX)], ejeY.tolist(), linewidth=width_linea, alpha=alfa_linea, color = linea_color, zorder=0)
                AX2.add_line(line)
                
                width_f = size_point
                factor1 = (tam_plot1*fig_per)/tam_plot1 # determina un primer factor a partir de las dimensiones del grafico x/y
                factor2 = (AX2.get_ylim()[1] - AX2.get_ylim()[0])/AX2.get_xlim()[1] # determina un factor a partir de los limites de los ejes 'X' y 'Y'
                factor3 = width_f*factor1*factor2 # producto de los factores mas el ancho o radio de las figuras
                
                centroY = (AX2.get_ylim()[1] - AX2.get_ylim()[0]) / 2
                
                ejexpos = []
                labelx = []
                for e, ejes in enumerate(zip(ejeX, ejeY, ColormaP)):
                    if (ejes[1] % 1) == 0:
                        valorY = int(ejes[1]) # numeros sin decimales son conviertidos a enteros
                    else:
                        valorY = round(ejes[1], 3)
                    
                    if markers_point == '⏺':
                        AX2.add_patch(Ellipse((e, valorY), width_f, factor3, angle= 0, linewidth=0, facecolor= ejes[2], zorder=2, alpha = alfa_marker,
                                                                  label = ejes[0], url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[ejes[0]]))
                    if markers_point == '⏹':
                        AX2.add_patch(Rectangle((e - (width_f/2), valorY - (factor3/2)), width_f, factor3, angle= 0, linewidth=0, facecolor= ejes[2], alpha = alfa_marker,
                                                                    zorder=2, label = ejes[0], url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[ejes[0]]))

                    if valorY < centroY:
                        AX2.text(e+0.08, valorY, '    '+str(round(valorY, 3)), color = text_color, ha = 'center', va = 'bottom', rotation = 90, fontfamily = family1, fontsize = text_size, alpha = text_alfa)
                    if valorY > centroY:
                        AX2.text(e+0.08, valorY, str(round(valorY, 3))+'    ', color = text_color, ha = 'center', va = 'top', rotation = 90, fontfamily = family1, fontsize = text_size, alpha = text_alfa)


                    ejexpos.append(e)
                    labelx.append(ejes[0])
                
                AX2.set_xticks(ejexpos)
                AX2.set_xticklabels(labelx, rotation = 90, ha= 'center', va = 'top', fontsize = ticklabels_X, fontname=family_axis)
                
                
                AX2.set_xlabel('Sample', fontsize=label_X, fontname=family_axis, weight="bold")
                
                AX2.set_ylabel(tipo_indice, fontsize=label_Y, fontname=family_axis, weight="bold")
                

                for labejey in AX2.yaxis.get_ticklabels():
                    labejey.set_fontsize(ticklabels_Y)

                 
                for ticky in AX2.yaxis.get_major_ticks():
                    ticky.label.set_fontfamily([family_axis])

                
                
                AX2.set_xticklabels(KITS[tipo_kits_1])
                
                if cambiar_sample == 'Code1':
                    etiquetas = KITS[tipo_kits_1]
                    AX2.set_xticklabels(etiquetas)
                if cambiar_sample == 'Code2':
                    etiquetas = [name_code[i] for i in KITS[tipo_kits_1]]
                    AX2.set_xticklabels(etiquetas)
                if cambiar_sample == 'Code3':
                    etiquetas = [name_code2[i] for i in KITS[tipo_kits_1]]
                    AX2.set_xticklabels(etiquetas)
                
                
                display(fig2)
                
            set_locus()
            
            
        out_updata_plot = widgets.interactive_output(updata_plot, {'tipo_kits_1':tipo_kits_1, 'tipo_indice':tipo_indice, 'tam_plot1':tam_plot1,
                                                                  'width_linea':width_linea, 'text_size':text_size, 'linea_color':linea_color, 'alfa_linea':alfa_linea,
                                                                   'size_point':size_point, 'markers_point':markers_point,
                                                                   'family1':family1, 'label_X':label_X, 'ticklabels_X':ticklabels_X,
                                                                   'label_Y':label_Y, 'ticklabels_Y':ticklabels_Y, 'family_axis':family_axis, 'text_color':text_color,
                                                                   'text_alfa':text_alfa, 'cambiar_sample':cambiar_sample, 'alfa_marker':alfa_marker, 'Multiple_colorS':Multiple_colorS,
                                                                   'ajustar_ejes':ajustar_ejes})
        
        
        
        
        items_plot_1 = VBox([tam_plot1,
                            width_linea,
                            alfa_linea,
                            text_size,
                            HBox([blanca, widgets.Label('Text font:'), family1]),
                            text_alfa,
                            size_point,
                             alfa_marker,
                            HBox([blanca, widgets.Label('Markers:'), markers_point, widgets.Label('Line color:'), linea_color, widgets.Label('Text color:'), text_color])
                            ])
        
        items_plot_1_box = Box(children=[items_plot_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_plot_1.children) * 31))+'px'))
        
        items_axis_1 = VBox([label_X, ticklabels_X, label_Y, ticklabels_Y, HBox([blanca, widgets.Label('Axis font:'), family_axis]), HBox([blanca, widgets.Label('Rename samples:'), cambiar_sample])])
        items_axis_1_box = Box(children=[items_axis_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_axis_1.children) * 31))+'px'))
        
        items_save_1 = VBox([HBox([widgets.HTML('<font color = grey> <b style="font-size:0.7vw">SAVE CHART: </b>'), blanca, filename_plot]),
                             HBox([blanca, widgets.Label('Formats:'), png1, jpeg1, svg1, pdf1]),
                     HBox([blanca, widgets.Label('Chart parameters:'), filename, params])])
        items_save_1_box = Box(children=[items_save_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_save_1.children) * 34))+'px'))
        
        
        
        items_data = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'DATA']
        items_axis = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'AXIS']
        items_save = [widgets.HTML('<font color = white> <b style="font-size:0.7vw">'+i+'</b>') for i in 'S']
        
        
        IndiceS = HBox([VBox([widgets.HTML('<font color = #1976d2> <b style="font-size:0.8vw">KITS</b>'), tipo_kit_1_box, #blanca,
                              widgets.HTML('<font color = grey> <b style="font-size:0.8vw">INDICES</b>'), tipo_indice_box,
                              VBox([widgets.HTML('<font color = gray> <h style="font-size:0.55vw">Adjust Y axis to min value:</h>'), ajustar_ejes])]), blanca,
                        
                       VBox([widgets.HTML('<font color = grey> <i class="fa fa-cog fa-2x fa-fw"></i> <b style="font-size:0.8vw">PLOT SETTINGS</b>'),
                             
                             HBox([Box(children =[VBox(items_data)]), items_plot_1_box]),
                             
                             HBox([Box(children =[VBox(items_axis)]), items_axis_1_box]),
                             
                             HBox([Box(children =[VBox(items_save)]), items_save_1_box])
                            ]),
                       VBox([widgets.HTML('<font color = grey> <b style="font-size:0.6vw">MARKER COLORS</b>'), Multiple_colorS, HBox([blanca, blanca, ayuda3])]), blanca, VBox([out_updata_plot])])

        display(IndiceS)

                
INDICES_button.on_click(button_clicked)
IND_DIV = VBox([INDICES_button, INDICES_output])





DIVERSITY_CALCULATIONS = VBox([PHY_DIV, IND_DIV])














# # Richness




import itertools





colores = {}
for q in qualitative_colors:
    colores[q] = [matplotlib.colors.to_hex(i) for i in plt.get_cmap(q)(np.arange(qualitative_colors[q]))]











"""
este diccionari contiene cada variable y dentro los elementos ordenados ya sea por orden alfabetico o numerico
"""
min_max_variables = []
variables_items_ordened = {}
for i in variables:
    df = metadata[i].unique()
    min_max_variables.append(len(df))
    ListA = df.tolist()
    #---
    try:
        ListA.sort(key=float)
    except ValueError:
        ListA.sort(key=str)
    variables_items_ordened[i] = ListA











"""
Este diccionario contiene las rampas de colores con posibilidad de generar una lista de longitud n, de acuerdo al numero de elementos dentro de cada variable.
"""
diccionarios_colores = {}
for m in [1] + sorted(set(min_max_variables)):
    colores2 = {}
    for c in colores:
        rampa = c
        num = 0
        num2 = m
        for e, i in enumerate(colores[rampa]):
            if len(colores[rampa][num: num2]) == m:
                #barcolor(lista = colores[rampa][num: num2])
                colores2[rampa+'_v'+str(e)] = colores[rampa][num: num2]
                num += 1
                num2 += 1
            else:
                pass
    diccionarios_colores[m] = colores2











from matplotlib.patches import Rectangle

def barcolor_v1(lista1 = [], lista2 = []):
    mpl.rcParams.update(mpl.rcParamsDefault)
    fig = plt.figure(frameon=False)

    #plt.style.use('seaborn') # plt.style.available
    #plt.rcParams['savefig.facecolor'] = 'b'

    largo = 1
    ancho = 0.2
    ax1 = fig.add_axes([0, 0.08, ancho, 0.1])
    
    
    rango2 = {}
    factor2 = largo/len(lista2)
    #print('#####', factor2)
    num = 0
    yy = factor2/2
    for i, j in zip(lista2, list(range(len(lista2)))):
        ax1.add_patch( Rectangle((num, 0), 
                                factor2, 0.3, 
                                fc =i,  
                                ec ='none', 
                                lw = 0))
        num += factor2
        rango2[i] = [j, yy]
        yy += factor2
        
        
    #num2 = factor2/2
    #for r in rango2.values():
    #    ax1.text(num2, 0.3, '', ha = 'center', va = 'center', fontsize = 8)
    #    num2 += factor2
    #ax1.text(0, 0.15, 'Variant ', ha = 'right', va = 'center', fontsize = 9)
    
    
    
    #####
    
    
    rango1 = {}
    factor1 = largo/len(lista1)
    #print('#####', factor1)
    num = 0
    xx = factor1/2
    for i, j in zip(lista1, list(range(len(lista1)))):
        ax1.add_patch( Rectangle((num, 0.7), 
                                factor1, 0.3, 
                                fc =i,  
                                ec ='none', 
                                lw = 0))
        num += factor1
        rango1[i] = [j, xx]
        xx += factor1
        
        
    #num2 = factor1/2
    #for r in rango1.values():
    #    ax1.text(num2, 1, '', ha = 'center', va = 'center', fontsize = 8)
    #    num2 += factor1
        
    #ax1.text(0, 0.85, 'Reference ', ha = 'right', va = 'center', fontsize = 9)


    
    
    for k in rango2:
        #print(rango2[k], rango1[k])
        #             x1            x2           y1    y2
        ax1.plot([rango2[k][1], rango2[k][1]], [0.29, 0.29+0.1], color = k, linestyle='solid', lw=1)
        ax1.plot([rango1[k][1], rango1[k][1]], [0.71, 0.71-0.1], color = k, linestyle='solid', lw=1)
        
        ax1.plot([rango2[k][1], rango1[k][1]], [0.29+0.1, 0.71-0.1], color = k, linestyle='solid', lw=1)
        
    
    
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])
    ax1.set_xlim(0, largo)
    ax1.set_ylim(0, largo)
    ax1.axis('off')


    plt.show()
    
    
def barcolor_v2(lista1 = []):
    mpl.rcParams.update(mpl.rcParamsDefault)
    fig = plt.figure(frameon=False, facecolor = 'red')

    #plt.style.use('seaborn') # plt.style.available
    #plt.rcParams['savefig.facecolor'] = 'b'

    largo = 1
    ancho = 0.20
    ax1 = fig.add_axes([0, 0, ancho, 0.04])
    
    rango1 = {}
    factor1 = largo/len(lista1)
    #print('#####', factor1)
    num = 0
    xx = factor1/2
    for i, j in zip(lista1, list(range(len(lista1)))):
        ax1.add_patch( Rectangle((num, 0), 
                                factor1, 1, 
                                fc =i,  
                                ec ='none', 
                                lw = 0))
        num += factor1
        rango1[i] = [j, xx]
        xx += factor1
        
    ax1.tick_params(bottom=False, right=False, top=False, left=True, width = 0, length=0, color='none')
    
    plt.xticks([])
    plt.yticks([]) 

    
    ax1.axis('off')


    plt.show()
    
def BarcolorR(color = ''):
    mpl.rcParams.update(mpl.rcParamsDefault)
    if color in list(qualitative_colors.keys()):
        palette = sns.color_palette(color, qualitative_colors[color])
        sns.palplot(palette)
        plt.gcf().set_size_inches(2.5, 0.27)
        plt.gca().axis('off')
        plt.show()
    else:
        palette = sns.color_palette(color, 30)
        sns.palplot(palette)
        plt.gcf().set_size_inches(2.5, 0.27)
        plt.gca().axis('off')
        plt.show()























blanca2 = widgets.Button(layout=Layout(width='10px', height='12px'), disabled=True)
blanca2.style.button_color = 'white'





figu = plt.figure(figsize=(3, 3))
ax = figu.add_axes([0, 0, 1, 1])
marcas = {}
for m in ['o', 's']: # for m in marker_dict:
    b = ax.scatter([1,2],[3,4], marker=m)
    square_mk, = b.get_paths()
    marcas[m] = square_mk
figu.clear()
plt.close()











RICHNESS_button = widgets.Button(description="PROCESS AND VISUALIZE", icon = 'fa-eye', layout=Layout(width='590px'))
RICHNESS_button.style.button_color = 'gainsboro' #'deepskyblue'
RICHNESS_button.style.font_weight = 'bold'
RICHNESS_output = widgets.Output()

family = sorted(['Liberation Serif','Microsoft Sans Serif','Open Sans','Times New Roman','3ds Light','Calibri','Comic Sans MS',
                  'Arial','Courier New','Microsoft Yi Baiti','Lucida Console'])

def button_clicked(b):
    
    import os.path
    import seaborn as sns
    
    if os.path.isfile('tablas/ASVs_PD.txt'):
        PD  = pd.read_csv('tablas/ASVs_PD.txt', sep = '\t')
        Indices = pd.read_csv('tablas/ASVs_Indices.txt', sep = '\t')
        Indices = pd.concat([Indices, PD]).reset_index(drop = True)
    else:
        Indices = pd.read_csv('tablas/ASVs_Indices.txt', sep = '\t')
    
    Indices_transpose = Indices.T.iloc[1:,:].astype(float)
    Indices_transpose.columns = list(Indices['Index'])
    Indices_transpose.insert(loc = 0, column='Name Sample', value = Indices.iloc[:,1:].columns)
    Indices_transpose = Indices_transpose.reset_index(drop = True)
    Indices_transpose = Indices_transpose.merge(metadata[['Name Sample'] + variables], on = 'Name Sample', how = 'left')
    
    with RICHNESS_output:
        clear_output(True)
        
        blanca3 = widgets.Button(layout=Layout(width='10px', height='10px'), disabled=True)
        blanca3.style.button_color = 'white'
        
        IndicE = widgets.ToggleButtons(options= list(reversed(list(dict(OrderedDict(Counter({i : len(i) for i in Indices.Index.tolist()}).most_common())).keys()))), value = 'Chao-1', button_style = '')
        IndicE.style.button_width = '170px'
        tipo_indice_box = Box(children=[VBox([IndicE])], layout=Layout(border='1px solid gainsboro', width='180px', height='380px'))
        
        tipo_kits_1 = widgets.ToggleButtons(options= ['Both kits'] + metadata[VARIABLE_KIT].unique().tolist(), value = 'Both kits', button_style = 'primary')
        tipo_kits_1.style.button_width = '170px'
        tipo_kits_1.style.font_weight = 'bold'
        tipo_kit_1_box = Box(children=[VBox([tipo_kits_1])], layout=Layout(border='1px solid #1976d2', width='180px', height='95px'))
        
        #if 'Kit' in variables:
            #variables.remove('Kit')
        VariablE = widgets.ToggleButtons(options= variables, value = variables[0], button_style = 'warning')
        VariablE.style.button_width = '160px'
        tipo_variable_box = Box(children=[VBox([VariablE])], layout=Layout(border='1px solid #ff9800', width='170px', height='208px'))
        
        
        SignificanciA = widgets.SelectionSlider(options=[0.01, 0.02, 0.03, 0.04, 0.05],value=0.05,disabled=False, description = 'P-value <:', layout=Layout(width='280px', height='25px'),
                                                continuous_update=False,orientation='horizontal',readout=True)
        mostrat_SignificanciA = widgets.ToggleButtons(options=['True', 'False'], value = 'True')
        mostrat_SignificanciA.style.button_width = '55px'
        
        n = 2.5
        cor_size2 = {}
        for e, i in enumerate(range(12)):
            cor_size2[e+1] = round(n, 2)
            n += 0.5
        
        ancho = widgets.SelectionSlider(options=list(cor_size2.keys()),value=5,disabled=False,
                                              description = 'Chart width:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        alto = widgets.SelectionSlider(options=list(cor_size2.keys())[:8],value=3,disabled=False,
                                              description = 'Chart heigth:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        alfa_box = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=1,disabled=False,
                                              description = 'Box alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        #markers_point = widgets.ToggleButtons(options=['o', 's', 'p', 'H', '*', 'D'], value = 'o')
        #markers_point.style.button_width = '26px'
        
        markers_point = widgets.ToggleButtons(options=['⏺', '⏹'], value = '⏺')
        markers_point.style.button_width = '31px'
        
        marker_color = widgets.ColorPicker(concise=False, value='#000000', disabled=False, layout = Layout(width='97px', height='25px'))

        marker_size = widgets.SelectionSlider(options=range(0, 101),value=50,disabled=False,
                                                      description = 'Marker size:',
                                                continuous_update=False,orientation='horizontal',readout=True,
                                                    layout=Layout(width='400px', height='25px'))

        marker_alfa = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=1,disabled=False,
                                                      description = 'Marker alpha:',
                                                continuous_update=False,orientation='horizontal',readout=True,
                                                   layout=Layout(width='400px', height='25px'))
        
        
        label_X = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'X label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_X = widgets.SelectionSlider(options=range(5, 31),value=10,disabled=False,
                                              description = 'X tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))

        label_Y = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'Y label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_Y = widgets.SelectionSlider(options=range(5, 31),value=10,disabled=False,
                                              description = 'Y tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        family_axis = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                   layout = Layout(width='290px', height='25px'))
        
        angulos = list(range(0, 361, 10))
        rotation = widgets.SelectionSlider(options=angulos,value=0,disabled=False,
                                              #description = 'X tick rotation:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='300px', height='25px'))
        
        width_box = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.5,disabled=False,
                                                      description = 'Box width:',
                                                continuous_update=False,orientation='horizontal',readout=True,
                                                   layout=Layout(width='400px', height='25px'))
        
        filename = widgets.Text(value='', placeholder='File name.txt', description='', disabled=False, layout = Layout(width='195px', height='25px'))
        
        filename_plot = widgets.Text(value='', placeholder='Chart name', description='', disabled=False, layout = Layout(width='270px', height='25px'))
        
        def updata_plot(tipo_kits_1, IndicE, VariablE, SignificanciA, mostrat_SignificanciA, width_box):
            
            if tipo_kits_1 == 'Both kits':
                estudio = Indices_transpose[['Name Sample', IndicE, VariablE]]
                pt = pd.pivot_table(estudio, values='Name Sample', index=VariablE, aggfunc=len).reset_index()
                pat = pt[pt['Name Sample'] > 1][VariablE].astype(str).tolist() # aquellas variables con valores menores o igual a 1 seran eliminadas
                estudio = estudio[estudio[VariablE].str.contains('^'+'$|^'.join(pat)+'$') == True]
                
            if tipo_kits_1 == 'DNPowerSoil':
                estudio = Indices_transpose[['Name Sample', IndicE, VariablE]]
                estudio = DataFrame(DNPowerSoil, columns = ['Name Sample']).merge(estudio, on = 'Name Sample', how = 'left').reset_index(drop = True)
                pt = pd.pivot_table(estudio, values='Name Sample', index=VariablE, aggfunc=len).reset_index()
                pat = pt[pt['Name Sample'] > 1][VariablE].astype(str).tolist() # aquellas variables con valores menores o igual a 1 seran eliminadas
                estudio = estudio[estudio[VariablE].str.contains('^'+'$|^'.join(pat)+'$') == True]
                
            if tipo_kits_1 == 'DNMicrobial':
                estudio = Indices_transpose[['Name Sample', IndicE, VariablE]]
                estudio = DataFrame(DNMicrobial, columns = ['Name Sample']).merge(estudio, on = 'Name Sample', how = 'left').reset_index(drop = True)
                pt = pd.pivot_table(estudio, values='Name Sample', index=VariablE, aggfunc=len).reset_index()
                pat = pt[pt['Name Sample'] > 1][VariablE].astype(str).tolist() # aquellas variables con valores menores o igual a 1 seran eliminadas
                estudio = estudio[estudio[VariablE].str.contains('^'+'$|^'.join(pat)+'$') == True]
                
                
            
            VariablE_orderned = estudio[VariablE].unique().tolist()
            #---
            try:
                VariablE_orderned.sort(key=float)
            except ValueError:
                VariablE_orderned.sort(key=str)
            #---
            orden_in_plot = dict(zip(VariablE_orderned, range(len(VariablE_orderned))))
            
            
            datos = {}
            for E in VariablE_orderned:
                datos[E] = estudio[estudio[VariablE] == E][IndicE].tolist()

            combinaciones = list(itertools.combinations(VariablE_orderned, 2))

            Significativos = {}
            for d in combinaciones:
                PVALUE = stats.ttest_ind(datos[d[0]], datos[d[1]], equal_var=False).pvalue
                if PVALUE < SignificanciA:
                    if PVALUE < 0.009999:
                        Significativos[d[0]+'_vs_'+d[1]] = {'ini_end':[orden_in_plot[d[0]], orden_in_plot[d[0]], orden_in_plot[d[1]], orden_in_plot[d[1]]],
                                                      'pval':format(PVALUE, '.3e')}
                    else:
                        Significativos[d[0]+'_vs_'+d[1]] = {'ini_end':[orden_in_plot[d[0]], orden_in_plot[d[0]], orden_in_plot[d[1]], orden_in_plot[d[1]]],
                                                      'pval':format(PVALUE, '.3e')}
            
            
            
            mpl.rcParams.update(mpl.rcParamsDefault)
            sns.set(style="white")

            fig = plt.figure(figsize=(3, 3))

            AX3 = fig.add_axes([0, 0, 1, 1])


            AX3 = sns.boxplot(x=VariablE, y=IndicE, data=estudio, hue=None, saturation=1, width=width_box,
                             order = VariablE_orderned, linewidth=1, fliersize = 0)

            AX3.tick_params(bottom=True, right=False, top=False, left=True, width = 1, length=4, color='black')

            AX3 = sns.swarmplot(x=VariablE, y=IndicE, data=estudio, color="black", order = VariablE_orderned)

            labels_Y = AX3.get_yticks()

            limy = AX3.get_ylim()[1] + (AX3.get_ylim()[1] * 0.05)
            
            AX3.xaxis.get_label().set_weight("bold")
            AX3.yaxis.get_label().set_weight("bold")
            
            plt.xlabel(VariablE, weight="bold")
            plt.ylabel(IndicE, weight="bold")
            
            if mostrat_SignificanciA == 'True':
                if len(Significativos) == 0:
                    pass
                else:
                    distancia = limy * 0.07 # 8%
                    for i in Significativos:
                        #print(Significativos[i]['ini_end'], [limy-(limy * 0.02), limy, limy, limy-(limy * 0.02)])
                        AX3.plot(Significativos[i]['ini_end'], [limy-(limy * 0.02), limy, limy, limy-(limy * 0.02)], '-', linewidth=1, color='black')
                        AX3.text(Significativos[i]['ini_end'][-1], limy, '  '+str(Significativos[i]['pval']), ha = 'left', va = 'center', fontsize = 10, fontfamily = 'Open Sans')
                        limy += distancia

                    AX3.set_yticks(labels_Y)
            else:
                pass

            plt.gca().spines['right'].set_color(None)
            plt.gca().spines['top'].set_color(None)


            plt.close()
            
            
            
            
            
            colores3 = widgets.Dropdown(options= list(diccionarios_colores[len(VariablE_orderned)].keys()), value='tab20_v0',
                             layout=Layout(width='100px', height='28px'))
            

            # --- save
            ancho_plot_save = str(71)
            png1 = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
            png1.style.button_color = 'gold'
            output1 = widgets.Output()
            def button_clicked1(b):
                with output1:
                    clear_output(True)
                    if filename_plot.value is '':
                        nombre_grafico = 'Richness_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                    if filename_plot.value is not '':
                        nombre_grafico = filename_plot.value
                        
                    fig.savefig('plots_asv/richness/'+nombre_grafico+'.png', dpi = 900, bbox_inches= 'tight')
                    
                    
            png1.on_click(button_clicked1)
            #----
            jpeg1 = widgets.Button(description="JPEG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
            jpeg1.style.button_color = 'gold'
            output2 = widgets.Output()
            def button_clicked2(b):
                with output2:
                    clear_output(True)
                    if filename_plot.value is '':
                        nombre_grafico = 'Richness_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                    if filename_plot.value is not '':
                        nombre_grafico = filename_plot.value
                        
                    fig.savefig('plots_asv/richness/'+nombre_grafico+'.jpeg', dpi = 900, bbox_inches= 'tight')
                    
            jpeg1.on_click(button_clicked2)
            #----
            svg1 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
            svg1.style.button_color = 'gold'
            output3 = widgets.Output()
            def button_clicked3(b):
                with output3:
                    clear_output(True)
                    if filename_plot.value is '':
                        nombre_grafico = 'Richness_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                    if filename_plot.value is not '':
                        nombre_grafico = filename_plot.value
                        
                    fig.savefig('plots_asv/richness/'+nombre_grafico+'.svg', dpi = 900, bbox_inches= 'tight')
                    
            svg1.on_click(button_clicked3)
            #----
            pdf1 = widgets.Button(description="PDF", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
            pdf1.style.button_color = 'gold'
            output4 = widgets.Output()
            def button_clicked4(b):
                with output4:
                    clear_output(True)
                    if filename_plot.value is '':
                        nombre_grafico = 'Richness_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                    if filename_plot.value is not '':
                        nombre_grafico = filename_plot.value
                        
                    fig.savefig('plots_asv/richness/'+nombre_grafico+'.pdf', dpi = 900, bbox_inches= 'tight')
                    
            pdf1.on_click(button_clicked4)
            ### ---
            
            
            params = widgets.Button(description="SAVE", icon = 'fa-save', layout=Layout(width=ancho_plot_save+'px'))
            params.style.button_color = 'cyan'
            output5 = widgets.Output()
            def button_clicked5(b):
                with output5:
                    clear_output(True)
                    
                    if filename.value is '':
                        nombre_archivo = 'Richness_chart_params_ASVs.txt'
                    if filename.value is not '':
                        nombre_archivo = filename.value
                    
                    with open('plots_asv/richness/'+nombre_archivo, 'w') as fq:
                        fq.write('#Saved parameters\n')
                        fq.write('Clustering:ASVs\n' )
                        """
                        diccionario para guardar en un archivo los parametros usados para volver a reproducirlos posteriormente
                        """
                        parametros = {'Kit:':tipo_kits_1,
                                      'Index:':IndicE,
                                      'Variable:':VariablE,
                                      'Significance:':SignificanciA,
                                      'Show sig:':mostrat_SignificanciA,
                                      'Box colors:':colores3.value,
                                      'Chart width:':ancho.value,
                                      'Chart heigth:':alto.value,
                                      'Box width:':width_box,
                                      'Box alpha:':alfa_box.value,
                                      'Marker color:':marker_color.value+' | '+'('+', '.join([r+':'+str(int(255*i)) for i, r in zip(to_rgb(marker_color.value), ['RED', 'GREEN', 'BLUE'])])+')',
                                      'Marker:':{'⏺':'circle', '⏹':'square'}[markers_point.value],
                                      'Marker size:':marker_size.value,
                                      'Marked alpha:':marker_alfa.value,
                                      'X label:':label_X.value,
                                      'X tick label:':ticklabels_X.value,
                                      'Y label:':label_Y.value,
                                      'Y tick labels:':ticklabels_Y.value,
                                      'Axis font:':family_axis.value,
                                      'X tick rotation:':rotation.value}
                        for w in parametros:
                            fq.write(w+str(parametros[w])+'\n')

            params.on_click(button_clicked5)
            #------
            
            
            def show_rampa(colores3):
                patron = colores3.split('_v')[0]
                barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(VariablE_orderned)][colores3])
            OUTshow_rampa = widgets.interactive_output(show_rampa, {'colores3':colores3})
            
            
            def boxrampas3(colores3, ancho, alto, alfa_box, markers_point, marker_color, marker_size, marker_alfa,
                           label_X, ticklabels_X, label_Y, ticklabels_Y, family_axis, rotation):

                def set_locus():
                    
                    fig.set_size_inches((cor_size2[ancho], cor_size2[alto]))
                    
                    [i.set_facecolor(c) for i, c in zip(AX3.artists, diccionarios_colores[len(VariablE_orderned)][colores3])]
                        
                    #limy = AX3.get_ylim()[1] + (AX3.get_ylim()[1] * 0.05)
                    
                    #labels_Y = AX3.get_yticks()

                    # configuracion del scatter
                    for i in AX3.collections:
                        i.set_paths([marcas[{'⏺':'o', '⏹':'s'}[markers_point]]]) # forma de puntos
                    for i in AX3.collections:
                        i.set_alpha(marker_alfa) # alfa de puntos
                    for i in AX3.collections:
                        i.set_sizes([marker_size]) # size de puntos
                    
                    for i in AX3.collections:
                        i.set_color(marker_color) # alfa de puntos
                        
                        
                        

                    [i.set_alpha(alfa_box) for i in AX3.artists]

                    #[i.set_alpha(alfa_box) for i in AX3.artists]
                    
                    
                    # EJES
                    
                    [i.set_fontsize(ticklabels_X) for i in AX3.xaxis.get_ticklabels()]
                    [i.set_fontsize(ticklabels_Y) for i in AX3.yaxis.get_ticklabels()]
                    
                    [tickx.label.set_fontfamily([family_axis]) for tickx in AX3.xaxis.get_major_ticks()]
                    [tickx.label.set_fontfamily([family_axis]) for tickx in AX3.yaxis.get_major_ticks()]
                    
                    AX3.xaxis.get_label().set_fontsize(label_X)
                    AX3.yaxis.get_label().set_fontsize(label_Y)
                    AX3.xaxis.get_label().set_fontfamily([family_axis])
                    AX3.yaxis.get_label().set_fontfamily([family_axis])
                    
                    
                    
                    [i.set_rotation(rotation) for i in AX3.get_xticklabels()]

                    display(fig)
                    
                        

                set_locus()
            
            OUTboxrampas3 = widgets.interactive_output(boxrampas3, {'colores3':colores3, 'ancho':ancho, 'alto':alto, 'alfa_box':alfa_box,
                                                                    'markers_point':markers_point, 'marker_color':marker_color, 'marker_size':marker_size,
                                                                     'marker_alfa':marker_alfa, 'label_X':label_X, 'ticklabels_X':ticklabels_X, 'label_Y':label_Y,
                                                                     'ticklabels_Y':ticklabels_Y, 'family_axis':family_axis, 'rotation':rotation})            
            
            
            #------------------------------------
            items_save_1 = VBox([HBox([widgets.HTML('<font color = grey> <b style="font-size:0.7vw">SAVE CHART: </b>'), blanca, filename_plot]),
                                 HBox([blanca, widgets.Label('Formats:'), png1, jpeg1, svg1, pdf1]),
                         HBox([blanca, widgets.Label('Chart parameters:'), filename, params])])
            items_save_1_box = Box(children=[items_save_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_save_1.children) * 34))+'px'))
            
            
            display(VBox([HBox([ VBox([blanca2, widgets.HTML('<font color = grey> <b style="font-size:0.7vw">BOX COLORS</b>'), colores3]), VBox([blanca3, OUTshow_rampa]), VBox([items_save_1_box])]), OUTboxrampas3, output4]))
                
                        
        out_updata_plot = widgets.interactive_output(updata_plot, {'tipo_kits_1':tipo_kits_1, 'IndicE':IndicE, 'VariablE':VariablE,
                                                                  'SignificanciA':SignificanciA, 'mostrat_SignificanciA':mostrat_SignificanciA, 'width_box':width_box}) 
        
        items_plot_1 = VBox([ancho,
                             alto,
                             width_box,
                             alfa_box,
                             HBox([blanca, widgets.Label('Marker color:'), marker_color, widgets.Label('Marker:'), markers_point]),
                             marker_size,
                             marker_alfa])
        
        items_plot_1_box = Box(children=[items_plot_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_plot_1.children) * 31))+'px'))
        
        significancia_box = Box(children=[HBox([SignificanciA, mostrat_SignificanciA])], layout=Layout(border='1px solid gainsboro', width='410px', height='38px'))
        
        items_axis_1 = VBox([label_X, ticklabels_X, label_Y, ticklabels_Y, HBox([blanca, widgets.Label('Axis font:'), family_axis]), HBox([widgets.Label('X tick rotation:'), rotation])])
        items_axis_1_box = Box(children=[items_axis_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_axis_1.children) * 31))+'px'))
        
        
        
        items_data = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'DATA']
        items_axis = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'AXIS']
        items_save = [widgets.HTML('<font color = white> <b style="font-size:0.7vw">'+i+'</b>') for i in 'S']
        
        
        RichnesS = HBox([VBox([widgets.HTML('<font color = #1976d2> <b style="font-size:0.8vw">KITS</b>'), tipo_kit_1_box,
                              HBox([VBox([widgets.HTML('<font color = grey> <b style="font-size:0.8vw">INDICES</b>'), tipo_indice_box]),
                                    VBox([widgets.HTML('<font color = grey> <b style="font-size:0.8vw">VARIABLES</b>'), tipo_variable_box])
                                   ])]),blanca,
                         VBox([widgets.HTML('<font color = grey> <i class="fa fa-cog fa-2x fa-fw"></i> <b style="font-size:0.8vw">PLOT SETTINGS</b>'),
                               
                               HBox([Box(children =[VBox(items_save)]), significancia_box]),
                               
                              HBox([Box(children =[VBox(items_data)]), items_plot_1_box]),
                             
                              HBox([Box(children =[VBox(items_axis)]), items_axis_1_box])
                              
                              ]),
                       blanca, VBox([out_updata_plot])])
        
        
        display(RichnesS)
        
RICHNESS_button.on_click(button_clicked)





RICHNESS_CALCULATIONS = VBox([RICHNESS_button, RICHNESS_output])














# # TAXONOMY










def _remove_dups(L):
    """
    Remove duplicates AND preserve the original order of the elements.
    The set class is not guaranteed to do this.
    """
    seen_before = set([])
    L2 = []
    for i in L:
        if i not in seen_before:
            seen_before.add(i)
            L2.append(i)
    return L2





import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import *


def SELECTSAM(SAM_SELECT = ''):
    sumary111 = []
    with open('plots_asv/taxonomy/AAAbundancEEE.txt', 'r') as fq:
        for enu, line in enumerate(fq):
            line = line.rstrip()
            if enu == 0:
                header = line.split('\t')
            else:
                sumary111.append(line.split('\t'))
    AAAbundancEEE = DataFrame(sumary111, columns = ['Sample'] + header[1:])
    AAAbundancEEE = AAAbundancEEE.set_index('Sample')
    AAAbundancEEE = AAAbundancEEE.astype('float64')
    
    with open('plots_asv/taxonomy/dict_para_salvar.json', 'r') as fp:
        CorrepondenciA = json.load(fp)
        
    df_one = AAAbundancEEE.T[[SAM_SELECT]][AAAbundancEEE.T[[SAM_SELECT]][SAM_SELECT] > 0].sort_values(by =SAM_SELECT,ascending=False)
    
    df_two = metadata[metadata['Name Sample'] == SAM_SELECT]
    df_two = df_two.set_index('Name Sample').T
    
    namesmax = max([len(i) for i in df_one.index])+3
    nummax = max([len(str(round(i, 4))) for i in df_one[SAM_SELECT]])
    namesmax2 = max([len(i) for i in df_two.index])+2
    nummax2 = max([len(i) for i in df_two[SAM_SELECT]])
    
    Linaje = CorrepondenciA['Linaje']
    Data = CorrepondenciA['data']
    Porcentaje = CorrepondenciA['Porcentaje']
    
    import datetime
    import ctypes
    
    ctypes.windll.shcore.SetProcessDpiAwareness(1)

    root= tk.Tk()
    root.title("Sample Data")
    root.geometry("1100x470")
    root.configure(background='white')
    root.attributes("-topmost", True)


    cero = Label(root, text=SAM_SELECT, font=("Arial", 10,  "bold"), fg = 'white', bg = 'darkblue')
    cero.grid(column=0, row=0, sticky = W+E+S+N)


    cero1 = Label(root, text='16S Analysis', font=("Arial", 12,  "bold"), fg = 'silver', bg = 'white')
    cero1.grid(column=4, row=0, sticky = W+S)

    labebl = Label(root, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 0, row = 1)

    primera = Label(root, text= Linaje, font=("Arial", 8,  "bold"), fg = 'black', bg = 'gainsboro')
    primera.grid(column=0, row=2, sticky = W+E+S)

    segunda = Label(root, text='Relative abundance (%)', font=("Arial", 8,  "bold"), fg = 'black', bg = 'khaki')
    segunda.grid(column=1, row=2, sticky = W+E+S)


    #--------------------------------------------------------
    table = tk.Text(root, font=("Arial italic", 8), height=len(df_one), width=namesmax, fg = 'black', bg = 'antiquewhite')
    table.grid(column=0, row=3, sticky =  W+E+N)

    for h, i in enumerate(df_one.index):
        table.insert(tk.INSERT, '  '+i+' \n')
    table.config(state=DISABLED)

    table2 = tk.Text(root, font=("Arial", 8), height=len(df_one), width=nummax, fg = 'black', bg ='lightcyan' )
    table2.grid(column=1, row=3, sticky =  W+E+N)
    for h, i in enumerate(df_one[SAM_SELECT]):
        table2.insert(tk.INSERT, '  '+str(round(i, 4))+'\n')
    table2.config(state=DISABLED)


    ##########################################################
    labebl = Label(root, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 0, row = 4)


    primera = Label(root, text='Attribute ', font=("Arial", 8,  "bold"), fg = 'black', bg = 'gainsboro')
    primera.grid(column=0, row=5, sticky = W+E+S)

    segunda = Label(root, text='Information', font=("Arial", 8,  "bold"), fg = 'black', bg = 'khaki')
    segunda.grid(column=1, row=5, sticky = W+E+S)

    #--------------------------------------------------------
    table3 = tk.Text(root, font=("Arial", 8,  "bold"), height=len(df_two), width=namesmax2, fg = 'black', bg = 'antiquewhite')
    table3.grid(column=0, row=6, sticky =  W+E+N)

    for h, i in enumerate(df_two.index):
        table3.insert(tk.INSERT, '  '+i+'\n')
    table3.config(state=DISABLED)

    table4 = tk.Text(root, font=("Arial", 8), height=len(df_two), width=nummax2, fg = 'black', bg ='lightcyan' )
    table4.grid(column=1, row=6, sticky =  W+E+N)
    for h, i in enumerate(df_two[SAM_SELECT]):
        table4.insert(tk.INSERT, '  '+i+'\n')
    table4.config(state=DISABLED)



    #--------------------------------------------------------
    labebl = Label(root, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 3, row = 4)

    labebl2 = Label(root, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl2.grid(column = 0, row = 7)


    agrupar = LabelFrame(root, text = "Save chart:", font=("Arial", 8))
    agrupar.grid(column=4, row=0, sticky= W+E+N+S)
    agrupar.configure(background='white')

    labebl2 = Label(agrupar, text= 'Formats:', font=("Arial", 8), fg="gray", bg = 'white')
    labebl2.grid(column = 0, row = 0)

    def salvarpng():

        figure.savefig('plots_asv/taxonomy/'+SAM_SELECT+'_'+Data+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')



    boton1 = Button(agrupar, text=" PNG ", cursor="hand2", borderwidth=0,
                activebackground= 'gainsboro',
                bg="darkblue", fg="white",font=("Arial", 7),
                    command = salvarpng) # newwin.destroy
    boton1.grid(column = 1, row = 0)

    labebl = Label(agrupar, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 2, row = 0)

    def salvarjpeg():

        figure.savefig('plots_asv/taxonomy/'+SAM_SELECT+'_'+Data+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.jpeg', dpi = 900, bbox_inches= 'tight')



    boton2 = Button(agrupar, text=" JPEG ", cursor="hand2", borderwidth=0,
                activebackground= 'gainsboro',
                bg="darkblue", fg="white",font=("Arial", 7),
                    command = salvarjpeg) # newwin.destroy
    boton2.grid(column = 3, row = 0)

    labebl = Label(agrupar, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 4, row = 0)

    def salvarsvg():

        figure.savefig('plots_asv/taxonomy/'+SAM_SELECT+'_'+Data+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')

    boton3 = Button(agrupar, text=" SVG ", cursor="hand2", borderwidth=0,
                activebackground= 'gainsboro',
                bg="darkblue", fg="white",font=("Arial", 7),
                    command = salvarsvg) # newwin.destroy
    boton3.grid(column = 5, row = 0)

    labebl = Label(agrupar, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 6, row = 0)


    def salvarpdf():

        figure.savefig('plots_asv/taxonomy/'+SAM_SELECT+'_'+Data+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.pdf', dpi = 900, bbox_inches= 'tight')

    boton4 = Button(agrupar, text=" PDF ", cursor="hand2", borderwidth=0,
                activebackground= 'gainsboro',
                bg="darkblue", fg="white",font=("Arial", 7),
                    command = salvarpdf) # newwin.destroy
    boton4.grid(column = 7, row = 0)

    labebl2 = Label(agrupar, text= '    ', font=("Arial", 8), fg="red", bg = 'white')
    labebl2.grid(column = 8, row = 0)
    #-------------------------------------


    group_aspect = LabelFrame(root, text = "", font=("Arial", 8,  "bold"), height = 1)
    group_aspect.grid(column=4, row=2, rowspan = 100, sticky= E+N)
    group_aspect.configure(background='white')


    mpl.rcParams.update(mpl.rcParamsDefault)
    figure = plt.figure(figsize=(6,3))

    bar1 = FigureCanvasTkAgg(figure, group_aspect)
    bar1.draw()
    bar1.get_tk_widget().grid(row=0, column = 0)#, rowspan=7, columnspan = 7


    ax = figure.add_axes([0, 0, 1, 1])

    ax.pie(df_one[SAM_SELECT].tolist(), #autopct ='%1.1f%%',
           labels = None,
                                pctdistance = 1, labeldistance= 1,
                                startangle = 0, radius = 0.4, rotatelabels = True,frame=True,
           center=(0.5, 0.5),
                                colors = [CorrepondenciA[i] for i in df_one.index],
                                wedgeprops={'alpha':1, 'linewidth': 0, 'edgecolor':'black'},
                                explode = np.array([0.0]*len(df_one)), textprops=dict(size = 8))
    if len(df_one) <= 17:
        plt.legend(df_one.index, bbox_to_anchor= [0.47,0.9], loc = 2,
                   handletextpad=0.5,ncol= 1, title=Linaje, title_fontsize = 7,
                                   fancybox=True, framealpha=0.5, shadow=False,
                                   handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                                   borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                                   prop={'style':'italic', 'size':6.5})
    else:
        plt.legend(df_one.index, bbox_to_anchor= [0.47,0.9], loc = 2,
                   title=Linaje, title_fontsize = 7,
                   handletextpad=0.5,ncol= 2,
                                   fancybox=True, framealpha=0.5, shadow=False,
                                   handlelength = .7, labelspacing = 0.5, columnspacing = 1,
                                   borderpad = 0.5, edgecolor="gainsboro", #frameon=False,
                                   prop={'style':'italic', 'size':6.5})

    centre_circle = plt.Circle((0.5,0.5),0.2,fc = 'white')
    plt.gca().add_artist(centre_circle)

    ax.text(0.5,0.5, SAM_SELECT, fontsize = 8.5, weight='bold', color = 'black', ha = 'center', va = 'center', zorder = 2)

    ax.set_xlim(0,2)

    ax.axis('off')
    plt.close()

    root.mainloop()





from collections import Counter
import matplotlib
from matplotlib import cm
import matplotlib.path as mpath
import matplotlib.patches as mpatches


from tkinter import ttk



def source1(radio = 1, theta1 = 0, theta2 = 45, width = 0.05, sep = 0.05, color = 'red'):
    sour = mpatches.Wedge((0, 0), # centro
                          radio, # radio
                          theta1, # angulo inicial
                          theta2, # angulo final
                          width=width, # ancho el 10% del radio
                          #url = 'https://www.uniprot.org/',
                          facecolor=color, zorder = 1)
    mitad = int(len(sour.get_path().vertices[:-2])/2)
    mitad1 = sour.get_path().vertices[:-2][:mitad,:]
    mitad2 = sour.get_path().vertices[:-2][mitad:,:]
    mitad3 = np.array(list(reversed(mitad2)))
    
    centro1 = int((len(mitad1)/2) - 0.5)
    centro2 = int((len(mitad2)/2) - 0.5)
    internopos = mitad2[centro2]
    etiquetapos = mitad1[centro1]
    
    return sour, mitad3, etiquetapos, mitad1
def source2(radio = 1, theta1 = 0, theta2 = 45, width = 10, sep = 0.05, color = 'red'):

    #............................................
    sour2 = mpatches.Wedge((0, 0), # centro
                           radio*(1-(width/100)), # radio
                           theta1, # angulo inicial
                           theta2, # angulo final
                           width=sep, # ancho el 10% del radio
                           facecolor=color, zorder = 1, lw = 0, ec = color)
    mitadd = int(len(sour2.get_path().vertices[:-2])/2)
    mitad11 = sour2.get_path().vertices[:-2][:mitadd,:]
    mitad22 = sour2.get_path().vertices[:-2][mitadd:,:]
    mitad33 = np.array(list(reversed(mitad22)))
    
    centro11 = int((len(mitad11)/2) - 0.5)
    centro22 = int((len(mitad22)/2) - 0.5)
    internopos = mitad22[centro22]
    
    return sour2, mitad33, internopos
def LOCATIONS(source, target):
    """
    source: array con posiciones del origen
    target: arreglo con posiciones hacia el destino
    """
    Path = mpath.Path
    path_data = []
    path_data.append(tuple([Path.MOVETO, tuple(source[0])]))
    for i in source[1:-1]:
        path_data.append(tuple([Path.CURVE3, tuple(i)]))
    path_data.append(tuple([Path.LINETO, tuple(source[-1])]))
    path_data.append(tuple([Path.CURVE3, (0,  0)]))
    path_data.append(tuple([Path.LINETO, tuple(target[0])]))
    for j in target[1:-1]:
        path_data.append(tuple([Path.CURVE3, tuple(j)]))
    path_data.append(tuple([Path.LINETO, tuple(target[-1])]))
    path_data.append(tuple([Path.CURVE3, (0,  0)]))
    path_data.append(tuple([Path.LINETO, tuple(source[0])]))
    return path_data













def chord_plot():
    
    sumary111 = []
    uno = open('plots_asv/taxonomy/AAAbundancEEE.txt', 'r')
    for enu, line in enumerate(uno):
        line = line.rstrip()
        if enu == 0:
            header = line.split('\t')
        else:
            sumary111.append(line.split('\t'))
    uno.close()
    sumary2 = DataFrame(sumary111, columns = ['Name Sample'] + header[1:])
    sumary2 = sumary2.set_index('Name Sample')
    sumary2 = sumary2.astype('float64')

    with open('plots_asv/taxonomy/dict_para_salvar.json', 'r') as fp:
        CorrepondenciA = json.load(fp)

    title_kit = CorrepondenciA['data']

    with open('plots_asv/taxonomy/dict_variable_element_colors.json', 'r') as fp:
        dict_variable_element_colors = json.load(fp)

    sssuuu = sumary2
    if 'Others' in sssuuu.columns.tolist():
        sssuuu = sssuuu.drop(columns = 'Others')
    else:
        pass

    origenes = sssuuu.columns.tolist()
    destinos = sssuuu.index.tolist()
    net = []
    for j in sssuuu.columns:
        a = sssuuu[[j]]
        a = a[a[j] > 0].index.tolist()
        for x in a:
            net.append(tuple([j, x]))
    Counts_SOURCE = dict(Counter([SOURCE for SOURCE, TARGET in net]))
    #print(Counts_SOURCE)
    Counts_TARGET = dict(Counter([TARGET for SOURCE, TARGET in net]))
    #print(Counts_TARGET)
    pal_max_origenes = max([len(i) for i in origenes])
    pal_max_destinos = max([len(i) for i in destinos])

    samples_variables = DataFrame(sumary2.index.tolist(), columns = ['Name Sample']).merge(metadata, on = 'Name Sample', how = 'left')   

    Linaje = CorrepondenciA['Linaje']
    Data = CorrepondenciA['data']
    Porcentaje = CorrepondenciA['Porcentaje']


    import datetime
    import ctypes

    ctypes.windll.shcore.SetProcessDpiAwareness(1)


    root= tk.Tk()
    root.title("Sample Data")
    root.geometry("1000x750")
    root.configure(background='white')
    root.attributes("-topmost", True)


    cero = Label(root, text='Settings', font=("Arial", 10,  "bold"), fg = 'white', bg = 'darkblue')
    cero.grid(column=0, row=0, sticky = 'WES')

    SIZE_TEXT_ALL = 6

    text_font = ("Arial", SIZE_TEXT_ALL)

    root.option_add('*TCombobox*Listbox.font', text_font)

    #-------------------------------------
    def on_select(event):

        mpl.rcParams.update(mpl.rcParamsDefault)

        fig = plt.figure(figsize=(7, 7))

        bar1 = FigureCanvasTkAgg(fig, group_aspect)
        bar1.draw()
        bar1.get_tk_widget().grid(row=0, column = 0)#, rowspan=7, columnspan = 7


        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_aspect('equal', 'box')

        sources_label = 1.36
        targets_label = 1.58

        radio = float(radio_source.get())
        sepp = float(source_espacio.get())
        ancho = source_width.get()
        #tam = 5
        Espacio = float(source_interespacio.get())


        constante = opening_source.get() - len(origenes)
        tam = constante / len(origenes) # 

        W = tam * len(origenes)

        Q = (Espacio * len(origenes)) - Espacio
        sour_inicial = 360 - ((W + Q)/ 2)
        #sour_inicial = 360 - 90

        """
        Nodos de Sources sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        """
        nodos_ori_0 = {}
        dif0 = 360 - ((W + Q)/ 2)
        dif1 = dif0 # inicio de los nodos
        record_ori_0 = {}

        for o in origenes:
            teta1 = dif0
            teta2 = dif0+tam
            sour, mitad3, etiquetapos, mitad1 = source1(radio = radio, theta1 = teta1, theta2 = teta2,
                                                    width = radio*(ancho/100), sep = sepp, color = CorrepondenciA[o])
            ax.add_patch(sour)
            record_ori_0[o] = [teta1, teta2]

            central_angle0 = (((dif0+tam)-dif0)/2)+dif0
            central_angle = (central_angle0*np.pi)/180
            if etiquetapos[0] < 0:
                central_angle = central_angle - np.pi
                if len(o) < pal_max_origenes:
                    tam_pal = len(o)
                    palabra = ' '*(pal_max_origenes - tam_pal)+o
                    palabra = palabra+' '*len(palabra)+' '
                else:
                    palabra = o
                    palabra = palabra+' '*len(palabra)+' '
            else:
                central_angle = (central_angle0*np.pi)/180
                if len(o) < pal_max_origenes:
                    tam_pal = len(o)
                    palabra = o+' '*(pal_max_origenes - tam_pal)
                    palabra = ' '*len(palabra)+' '+palabra
                else:
                    palabra = o
                    palabra = ' '*len(palabra)+' '+palabra


            ax.text(etiquetapos[0], etiquetapos[1], palabra, color = 'black',
                    va = 'center', ha = 'center', #fontweight = 'bold',
                    fontsize = source_letra.get(), rotation = np.rad2deg(central_angle),
                    family = 'monospace', style='italic')
            #ax.scatter(etiquetapos[0], etiquetapos[1], s = 50, c = 'black', zorder = 2)

            if Counts_SOURCE[o] == 1:
                t1 = dif1
                t2 = dif1+tam
                sour2, mitad33, internopos = source2(radio = radio, theta1 = t1, theta2 = t2,
                                                    width = ancho, sep = sepp, color = 'white')
                ax.add_patch(sour2)
                nodos_ori_0[o] = [sour2, mitad33, internopos]
                dif1 += (tam+Espacio)
            else:
                t1 = dif1
                sectores = tam/Counts_SOURCE[o]
                ss0 = 0
                ss1 = sectores
                intersecciones_ori = []
                for r in range(Counts_SOURCE[o]):
                    t1 = ss0+dif1
                    t2 = ss1+dif1
                    sour2, mitad33, internopos = source2(radio = radio, theta1 = t1, theta2 = t2,
                                                        width = ancho, sep = sepp, color = 'white')

                    ax.add_patch(sour2)
                    intersecciones_ori.append([sour2, mitad33, internopos])

                    ss1 += sectores
                    ss0 += sectores
                nodos_ori_0[o] = intersecciones_ori
                dif1 += (tam+Espacio)

            PER = 2
            central_angle = (((dif0+tam)-dif0)/2)+dif0

            dif0 += (tam+Espacio)

        """
        Nodos de Targets ttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
        """   
        #----------
        # angulo de separacion elegido



        angulo_espacio = SEP.get()  #+++++++++++++++++++++

        continuacion = dif0 - Espacio - 360
        Espacio2 = float(target_interespacio.get())  #++++++++++++++++++++

        if angulo_espacio > 0:
            dif00 = angulo_espacio + continuacion
            tam2 = ((sour_inicial-angulo_espacio)-(continuacion+angulo_espacio) - ((len(destinos) * Espacio2)-Espacio2))/len(destinos)
        else:
            tam2 = ((sour_inicial-continuacion) - ((len(destinos) * Espacio2)-Espacio2))/len(destinos)
            dif00 = continuacion

        nodos_des_0 = {}

        ANGULOS = []
        nodos_des_1 = {}
        #dif00 = continuacion
        dif11 = dif00
        #tam2 = ((sour_inicial-continuacion) - ((len(destinos) * Espacio2)-Espacio2))/len(destinos)

        radio2 = float(radio_target.get()) #++++++++++++++++++++++++++++++

        sepp2 = float(target_espacio.get()) #++++++++++++++++++++++++++++++


        ancho2 = target_width.get() #++++++++++++++++++++++++++++++


        var_radios = {list(dict_variable_element_colors.keys())[0]:radio2}

        increase = sepp2+(radio2*(ancho2/100))
        for i in list(dict_variable_element_colors.keys())[1:]:
            var_radios.update({i:radio2+increase})
            increase += sepp2+(radio2*(ancho2/100))

        MitadeS = []


        for o in destinos:
            teta11 = dif11
            teta22 = dif11+tam2



            if 'None' in list(dict_variable_element_colors.keys()):

                sour, mitad2, etiquetapos, mitad1 = source1(radio = radio2, theta1 = teta11, theta2 = teta22,
                                                         width = radio2*(ancho2/100), sep = sepp2, color = 'silver')
                ax.add_patch(sour)
                ETIQUETAS = etiquetapos

                MitadeS.append([mitad1, mitad2])



            else:
                #### Location
                if 'Location' in list(dict_variable_element_colors.keys()):
                    sour, mitad2, etiquetapos, mitad1 = source1(radio = var_radios['Location'], theta1 = teta11, theta2 = teta22,
                                                            width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Location'][correspondencia_sam_vars[o]['Location']])
                    ax.add_patch(sour)
                    ETIQUETAS = etiquetapos

                    if o == destinos[0]:
                        MITADES = [mitad1, mitad2]
                        MitadeS.append(MITADES)

                ##### ug OTA/kg
                if 'ug OTA/kg' in list(dict_variable_element_colors.keys()):
                    sour0, mitad22, etiquetapos0, mitad11 = source1(radio = var_radios['ug OTA/kg'], theta1 = teta11, theta2 = teta22,
                                                            width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['ug OTA/kg'][correspondencia_sam_vars[o]['ug OTA/kg']])
                    ax.add_patch(sour0)
                    ETIQUETAS = etiquetapos0

                    if o == destinos[0]:
                        MITADES = [mitad11, mitad22]
                        MitadeS.append(MITADES)


                ##### Cultivation
                if 'Cultivation' in list(dict_variable_element_colors.keys()):
                    sour00, mitad222, etiquetapos00, mitad111 = source1(radio = var_radios['Cultivation'], theta1 = teta11, theta2 = teta22,
                                                            width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Cultivation'][correspondencia_sam_vars[o]['Cultivation']])
                    ax.add_patch(sour00)
                    ETIQUETAS = etiquetapos00

                    if o == destinos[0]:
                        MITADES = [mitad111, mitad222]
                        MitadeS.append(MITADES)

                ##### Coffee Variety
                if 'Coffee Variety' in list(dict_variable_element_colors.keys()):
                    sour000, mitad2222, etiquetapos000, mitad1111 = source1(radio = var_radios['Coffee Variety'], theta1 = teta11, theta2 = teta22,
                                                            width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Coffee Variety'][correspondencia_sam_vars[o]['Coffee Variety']])
                    ax.add_patch(sour000)
                    ETIQUETAS = etiquetapos000

                    if o == destinos[0]:
                        MITADES = [mitad1111, mitad2222]
                        MitadeS.append(MITADES)

                ##### Genomic DNA kit
                if 'Genomic DNA kit' in list(dict_variable_element_colors.keys()):
                    sour0000, mitad22222, etiquetapos0000, mitad11111 = source1(radio = var_radios['Genomic DNA kit'], theta1 = teta11, theta2 = teta22,
                                                            width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Genomic DNA kit'][correspondencia_sam_vars[o]['Genomic DNA kit']])
                    ax.add_patch(sour0000)
                    ETIQUETAS = etiquetapos0000

                    if o == destinos[0]:
                        MITADES = [mitad11111, mitad22222]
                        MitadeS.append(MITADES)


                ##### Drying Time (Days)
                if 'Drying Time (Days)' in list(dict_variable_element_colors.keys()):
                    sour00000, mitad222222, etiquetapos00000, mitad111111 = source1(radio = var_radios['Drying Time (Days)'], theta1 = teta11, theta2 = teta22,
                                                            width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Drying Time (Days)'][correspondencia_sam_vars[o]['Drying Time (Days)']])
                    ax.add_patch(sour00000)
                    ETIQUETAS = etiquetapos00000

                    if o == destinos[0]:
                        MITADES = [mitad111111, mitad222222]
                        MitadeS.append(MITADES)

                ##### Postharvest Processing
                if 'Postharvest Processing' in list(dict_variable_element_colors.keys()):
                    sour000000, mitad2222222, etiquetapos000000, mitad1111111 = source1(radio = var_radios['Postharvest Processing'], theta1 = teta11, theta2 = teta22,
                                                            width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Postharvest Processing'][correspondencia_sam_vars[o]['Postharvest Processing']])
                    ax.add_patch(sour000000)
                    ETIQUETAS = etiquetapos000000

                    if o == destinos[0]:
                        MITADES = [mitad1111111, mitad2222222]
                        MitadeS.append(MITADES)

            #####
            central_angle0 = (((dif00+tam2)-dif00)/2)+dif00
            central_angle = (central_angle0*np.pi)/180


            if etiquetapos[0] < 0:
                central_angle = central_angle - np.pi
                if len(o) < pal_max_destinos:
                    tam_pal = len(o)
                    palabra = ' '*(pal_max_destinos - tam_pal)+o
                    palabra = palabra+' '+' '*len(palabra)
                else:
                    palabra = o
                    palabra = palabra+' '+' '*len(palabra)

            else:
                central_angle = (central_angle0*np.pi)/180
                if len(o) < pal_max_destinos:
                    tam_pal = len(o)
                    palabra = o+' '*(pal_max_destinos - tam_pal)
                    palabra = ' '+' '*len(palabra)+palabra
                else:
                    palabra = o
                    palabra = ' '+' '*len(palabra)+palabra


            ANGULOS.append(dif00)
            ax.text(ETIQUETAS[0], ETIQUETAS[1], palabra, color = 'black', va = 'center', ha = 'center', 
                    fontsize = target_letra.get(), rotation = np.rad2deg(central_angle), fontfamily = 'Liberation mono', fontweight='bold')

            #ax.scatter(etiquetapos33[0], etiquetapos33[1], s = 20, c = 'lime', zorder = 2)



            if Counts_TARGET[o] == 1:
                t1 = dif11
                t2 = dif11+tam2
                sour2, mitad33, internopos = source2(radio = radio2, theta1 = t1, theta2 = t2,
                                                    width = ancho2, sep = sepp2, color = 'white')
                ax.add_patch(sour2)
                nodos_des_1[o] = [sour2, mitad33, internopos]
                dif11 += (tam2+Espacio2)
            else:
                t1 = dif11
                sectores = tam2/Counts_TARGET[o]
                ss0 = 0
                ss1 = sectores
                intersecciones_des = []
                for r in range(Counts_TARGET[o]):
                    t1 = ss0+dif11
                    t2 = ss1+dif11
                    sour2, mitad33, internopos = source2(radio = radio2, theta1 = t1, theta2 = t2,
                                                        width = ancho2, sep = sepp2, color = 'white')
                    ax.add_patch(sour2)
                    intersecciones_des.append([sour2, mitad33, internopos])
                    ss1 += sectores
                    ss0 += sectores
                nodos_des_1[o] = intersecciones_des
                dif11 += (tam2+Espacio2)

            central_angle2 = (((dif00+tam2)-dif00)/2)+dif00
            radian2 = (central_angle2*np.pi)/180
            x2 = (radio2 * (1-((ancho2/2)/100))) *  np.cos(radian2)
            y2 = (radio2 * (1-((ancho2/2)/100))) *  np.sin(radian2)
            #ax.scatter(x2, y2, s = 10, c = 'black', zorder = 2)

            dif00 += (tam2+Espacio2)

        """
        Conexiones ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        """
        #---------------
        # parte 1
        XX = []
        for ori in origenes:
            xu = 0
            if Counts_SOURCE[ori] == 1:
                for SOURCE, TARGET in net:
                    if ori == SOURCE:
                        if Counts_TARGET[TARGET] == 1:  # target uno
                            #print('>>>', SOURCE, TARGET)
                            path_data = LOCATIONS(nodos_ori_0[SOURCE][1], nodos_des_1[TARGET][1])

                            codes, verts = zip(*path_data)
                            path = mpath.Path(verts, codes)
                            patch = mpatches.PathPatch(path, facecolor=CorrepondenciA[ori], alpha=1, lw = None, ec = 'none', zorder = 0)
                            ax.add_patch(patch)
                        else: # target con mas de uno
                            XX.append([SOURCE, 'NA', TARGET])
            else:
                for SOURCE, TARGET in net:
                    if ori == SOURCE:
                        #print(SOURCE, xu, TARGET)
                        if Counts_TARGET[TARGET] == 1:
                            #print(SOURCE, xu, TARGET, '----')
                            path_data = LOCATIONS(nodos_ori_0[SOURCE][xu][1], nodos_des_1[TARGET][1])
                            codes, verts = zip(*path_data)
                            path = mpath.Path(verts, codes)
                            patch = mpatches.PathPatch(path, facecolor=CorrepondenciA[ori], alpha=1, lw = None, ec = 'none', zorder = 0)
                            ax.add_patch(patch)
                        else: # target con mas de uno
                            XX.append([SOURCE, xu, TARGET])
                        xu += 1


        #---------------
        # Parte 2
        output = []
        for SOURCE, P, TARGET in XX:
            if TARGET not in output:
                output.append(TARGET)

        for s in output:
            n = 0
            for SOURCE, P, TARGET in XX:
                if s == TARGET:
                    if P == 'NA':
                        path_data = LOCATIONS(nodos_ori_0[SOURCE][1], nodos_des_1[TARGET][n][1])

                        codes, verts = zip(*path_data)
                        path = mpath.Path(verts, codes)
                        patch = mpatches.PathPatch(path, facecolor=CorrepondenciA[SOURCE], alpha=1, lw = None, ec = 'none', zorder = 0)
                        ax.add_patch(patch)
                    else:

                        path_data = LOCATIONS(nodos_ori_0[SOURCE][P][1], nodos_des_1[TARGET][n][1])

                        codes, verts = zip(*path_data)
                        path = mpath.Path(verts, codes)
                        patch = mpatches.PathPatch(path, facecolor=CorrepondenciA[SOURCE], alpha=1, lw = None, ec = 'none', zorder = 0)
                        ax.add_patch(patch)
                    n += 1
        #
        if mostrar_vari.get() == 'True':  #+++++++++++++++++++++++++++++
            if target_width.get() > 3:  #++++++++++++++++++
                for i, j in zip(MitadeS, list(dict_variable_element_colors.keys())):
                    ax.text(np.sum(i[0][0][0]+i[1][0][0])/2, np.sum(i[0][0][1]+i[1][0][1])/2,
                            ' '*len(j)+' '+j, color = 'black', ha = 'center',va = 'center', #weight = 'bold',
                            fontsize = 6, rotation = ((continuacion+angulo_espacio)-90)-2, fontfamily = 'Liberation mono')

        #
        if mostrar_vari.get() == 'False':  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            pass

        ##

        if 'None' in list(dict_variable_element_colors.keys()):
            pass
        else:

            if mostrar_leyenda.get() == 'True':  #++++++++++

                yy = 1.8
                xx = -2.4
                for var in list(dict_variable_element_colors.keys()):
                    ax.text(xx-0.09, yy, var, ha = 'left', va = 'center', fontsize = leyenda_letra.get(), fontfamily = 'Arial')
                    yy -= 0.07
                    for e, j in enumerate(dict_variable_element_colors[var]):
                        ax.scatter(xx-0.07, yy, s = leyenda_size.get(), c = dict_variable_element_colors[var][j], marker = 'o')
                        ax.text(xx-0.05, yy, '  '+j, ha = 'left', va = 'center', fontsize = leyenda_letra.get(), rotation = 0, fontfamily = '3ds Light')
                        #xx +=0.07
                        yy -= 0.07
                    yy -= 0.05
            if mostrar_leyenda.get() == 'False':  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                pass


        #ax.text(-0.35, 1.9, title_kit, ha = 'center', va = 'center', fontsize = 10, weight = 'bold')


        ax.set_xlim(-float(lim_xy.get())-0.5, float(lim_xy.get()))  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ax.set_ylim(-float(lim_xy.get()), float(lim_xy.get()))  #++++++++++++++++++
        ax.axis('off')

        plt.close()


        def on_select2():
            fig.savefig('plots_asv/taxonomy/Chord_'+title_kit+'_'+Data+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')
        ###############
        ###############


        boton = Button(agrupar, text="PNG", cursor="hand2", borderwidth=0,
                    bg="gold", fg="black",font=("Arial", SIZE_TEXT_ALL), command = on_select2)
        boton.grid(column = 0, row = 0, sticky = 'WES')
        #boton.bind('<Button-1>', on_select2)
        
        labebl = Label(agrupar, text= '', font=("Arial", 8), fg="red", bg = 'white')
        labebl.grid(column = 1, row = 0)

        def on_select3():
            fig.savefig('plots_asv/taxonomy/Chord_'+title_kit+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.jpeg', dpi = 900, bbox_inches= 'tight')
        ###############
        ###############


        boton = Button(agrupar, text="JPEG", cursor="hand2", borderwidth=0,
                    bg="gold", fg="black",font=("Arial", SIZE_TEXT_ALL), command = on_select3)
        boton.grid(column = 2, row = 0, sticky = 'WES')
        #boton.bind('<Button-1>', on_select2)
        
        labebl = Label(agrupar, text= '', font=("Arial", 8), fg="red", bg = 'white')
        labebl.grid(column = 3, row = 0)

        #-------------------------------------
        #-------------------------------------
        def on_select4():
            fig.savefig('plots_asv/taxonomy/Chord_'+title_kit+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')
        ###############
        ###############


        boton = Button(agrupar, text="SVG", cursor="hand2", borderwidth=0,
                    bg="gold", fg="black",font=("Arial", SIZE_TEXT_ALL), command = on_select4)
        boton.grid(column = 4, row = 0, sticky = 'WES')
        #boton.bind('<Button-1>', on_select3)
        
        labebl = Label(agrupar, text= '', font=("Arial", 8), fg="red", bg = 'white')
        labebl.grid(column = 5, row = 0)
        #-------------------------------------
        #-------------------------------------
        def on_select5():
            fig.savefig('plots_asv/taxonomy/Chord_'+title_kit+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.pdf', dpi = 900, bbox_inches= 'tight')
        ###############
        ###############


        boton = Button(agrupar, text="PDF", cursor="hand2", borderwidth=0,
                    bg="gold", fg="black",font=("Arial", SIZE_TEXT_ALL), command = on_select5)
        boton.grid(column = 6, row = 0, sticky = 'WES')
        #boton.bind('<Button-1>', on_select3)



    ###############
    color_botones = 'gainsboro'
    color_letra_boton = 'black'
    ancho_boton = 5
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero1 = Label(root, text='Separation', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero1.grid(column=0, row=1, sticky = 'WES')
    tamanos = list(range(5,101))
    SEP = IntVar()
    separacion = ttk.Combobox(root, textvariable = SEP, font=("Arial", SIZE_TEXT_ALL), style='W.TCombobox',
                         values = tamanos, width=ancho_boton)
    separacion.grid(column=1, row=1, sticky= 'SW')
    separacion.bind('<<ComboboxSelected>>', on_select)
    separacion.current(45)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero2 = Label(root, text='Margins', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero2.grid(column=0, row=2, sticky = 'WES')
    limites = ['1.5', '2', '2.5', '3', '3.5', '4']
    lim_xy = StringVar()
    margins = ttk.Combobox(root, textvariable = lim_xy, font=("Arial", SIZE_TEXT_ALL),style='W.TCombobox',
                         values = limites, width=ancho_boton)
    margins.grid(column=1, row=2, sticky= 'SW')
    margins.bind('<<ComboboxSelected>>', on_select)
    margins.current(1)


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero1 = Label(root, text='Source settings', font=("Arial", SIZE_TEXT_ALL, 'bold'), fg = 'black', bg = 'lime')
    cero1.grid(column=0, row=3, sticky = 'WES')

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero1 = Label(root, text='Source opening', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero1.grid(column=0, row=4, sticky = 'WES')
    apertura = list(range(20,221))
    opening_source = IntVar()
    abertura = ttk.Combobox(root, textvariable = opening_source, font=("Arial", SIZE_TEXT_ALL),
                         values = apertura, width=ancho_boton)
    abertura.grid(column=1, row=4, sticky= 'SW')
    abertura.bind('<<ComboboxSelected>>', on_select)
    abertura.current(160)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero00 = Label(root, text='Source radio', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero00.grid(column=0, row=5, sticky = 'WES')
    radio_fuente = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','1.5']
    radio_source = StringVar()
    source_radio = ttk.Combobox(root, textvariable = radio_source, font=("Arial", SIZE_TEXT_ALL),
                         values = radio_fuente, width=ancho_boton)
    source_radio.grid(column=1, row=5, sticky= 'SW')
    source_radio.bind('<<ComboboxSelected>>', on_select)
    source_radio.current(9)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero3 = Label(root, text='Source width', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero3.grid(column=0, row=6, sticky = 'WES')
    ancho_sources = list(range(0,51))
    source_width = IntVar()
    source_ancho = ttk.Combobox(root, textvariable = source_width, font=("Arial", SIZE_TEXT_ALL),
                         values = ancho_sources, width=ancho_boton)
    source_ancho.grid(column=1, row=6, sticky= 'SW')
    source_ancho.bind('<<ComboboxSelected>>', on_select)
    source_ancho.current(7)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero4 = Label(root, text='Source fontsize', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero4.grid(column=0, row=7, sticky = 'WES')
    letra_source = list(range(2,21))
    source_letra = IntVar()
    source_fontsize = ttk.Combobox(root, textvariable = source_letra, font=("Arial", SIZE_TEXT_ALL),
                         values = letra_source, width=ancho_boton)
    source_fontsize.grid(column=1, row=7, sticky= 'SW')
    source_fontsize.bind('<<ComboboxSelected>>', on_select)
    source_fontsize.current(6)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero4 = Label(root, text='Source space', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero4.grid(column=0, row=8, sticky = 'WES')
    espacio_source = ['0', '0.01', '0.025', '0.05', '0.075', '0.1', '0.11', '0.125', '0.15', '0.175', '0.2']
    source_espacio = StringVar()
    source_space = ttk.Combobox(root, textvariable = source_espacio, font=("Arial", SIZE_TEXT_ALL),
                         values = espacio_source, width=ancho_boton)
    source_space.grid(column=1, row=8, sticky= 'SW')
    source_space.bind('<<ComboboxSelected>>', on_select)
    source_space.current(1)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero4 = Label(root, text='Source interspace', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero4.grid(column=0, row=9, sticky = 'WES')
    interespacio_source = np.round(np.arange(0, 3.1, 0.1), 2).tolist()
    source_interespacio = StringVar()
    source_interspace = ttk.Combobox(root, textvariable = source_interespacio, font=("Arial", SIZE_TEXT_ALL),
                         values = interespacio_source, width=ancho_boton)
    source_interspace.grid(column=1, row=9, sticky= 'SW')
    source_interspace.bind('<<ComboboxSelected>>', on_select)
    source_interspace.current(5)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



    cero1 = Label(root, text='Target settings', font=("Arial", SIZE_TEXT_ALL, 'bold'), fg = 'black', bg = 'lime')
    cero1.grid(column=0, row=10, sticky = 'WES')
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero00 = Label(root, text='Target radio', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero00.grid(column=0, row=11, sticky = 'WES')
    radio_destino = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','1.5']
    radio_target = StringVar()
    target_radio = ttk.Combobox(root, textvariable = radio_target, font=("Arial", SIZE_TEXT_ALL),
                         values = radio_destino, width=ancho_boton)
    target_radio.grid(column=1, row=11, sticky= 'SW')
    target_radio.bind('<<ComboboxSelected>>', on_select)
    target_radio.current(9)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero5 = Label(root, text='Target width', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero5.grid(column=0, row=12, sticky = 'WES')
    ancho_targets = list(range(1,51))
    target_width = IntVar()
    target_ancho = ttk.Combobox(root, textvariable = target_width, font=("Arial", SIZE_TEXT_ALL),
                         values = ancho_targets, width=ancho_boton)
    target_ancho.grid(column=1, row=12, sticky= 'SW')
    target_ancho.bind('<<ComboboxSelected>>', on_select)
    target_ancho.current(5)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero6 = Label(root, text='Target fontsize', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero6.grid(column=0, row=13, sticky = 'WES')
    letra_target = list(range(2,21))
    target_letra = IntVar()
    target_fontsize = ttk.Combobox(root, textvariable = target_letra, font=("Arial", SIZE_TEXT_ALL),
                         values = letra_source, width=ancho_boton)
    target_fontsize.grid(column=1, row=13, sticky= 'SW')
    target_fontsize.bind('<<ComboboxSelected>>', on_select)
    target_fontsize.current(5)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero4 = Label(root, text='Target space', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero4.grid(column=0, row=14, sticky = 'WES')
    espacio_target = ['0', '0.01', '0.025', '0.05', '0.075', '0.1', '0.11', '0.125', '0.15', '0.175', '0.2']
    target_espacio = StringVar()
    target_space = ttk.Combobox(root, textvariable = target_espacio, font=("Arial", SIZE_TEXT_ALL),
                         values = espacio_target, width=ancho_boton)
    target_space.grid(column=1, row=14, sticky= 'SW')
    target_space.bind('<<ComboboxSelected>>', on_select)
    target_space.current(1)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero4 = Label(root, text='Target interspace', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero4.grid(column=0, row=15, sticky = 'WES')
    interespacio_target = np.round(np.arange(0, 3.1, 0.1), 2).tolist()
    target_interespacio = StringVar()
    target_interspace = ttk.Combobox(root, textvariable = target_interespacio, font=("Arial", SIZE_TEXT_ALL),
                         values = interespacio_target, width=ancho_boton)
    target_interspace.grid(column=1, row=15, sticky= 'SW')
    target_interspace.bind('<<ComboboxSelected>>', on_select)
    target_interspace.current(5)




    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    cero1 = Label(root, text='Legend settings', font=("Arial", SIZE_TEXT_ALL, 'bold'), fg = 'black', bg = 'lime')
    cero1.grid(column=0, row=16, sticky = 'WES')
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero00 = Label(root, text='Show legend', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero00.grid(column=0, row=17, sticky = 'WES')
    leyenda = ['True', 'False']
    mostrar_leyenda = StringVar()
    ley = ttk.Combobox(root, textvariable = mostrar_leyenda, font=("Arial", SIZE_TEXT_ALL),
                         values = leyenda, width=ancho_boton)
    ley.grid(column=1, row=17, sticky= 'SW')
    ley.bind('<<ComboboxSelected>>', on_select)
    ley.current(0)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero6 = Label(root, text='Legend fontsize', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero6.grid(column=0, row=18, sticky = 'WES')
    letra_leyenda = list(range(2,16))
    leyenda_letra = IntVar()
    leyenda_fontsize = ttk.Combobox(root, textvariable = leyenda_letra, font=("Arial", SIZE_TEXT_ALL),
                         values = letra_leyenda, width=ancho_boton)
    leyenda_fontsize.grid(column=1, row=18, sticky= 'SW')
    leyenda_fontsize.bind('<<ComboboxSelected>>', on_select)
    leyenda_fontsize.current(5)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero6 = Label(root, text='Marker size', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero6.grid(column=0, row=19, sticky = 'WES')
    size_leyenda = list(range(2,51))
    leyenda_size = IntVar()
    ley_size = ttk.Combobox(root, textvariable = leyenda_size, font=("Arial", SIZE_TEXT_ALL),
                         values = size_leyenda, width=ancho_boton)
    ley_size.grid(column=1, row=19, sticky= 'SW')
    ley_size.bind('<<ComboboxSelected>>', on_select)
    ley_size.current(8)

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero00 = Label(root, text='Variable name', font=("Arial", SIZE_TEXT_ALL), fg = color_letra_boton, bg = color_botones)
    cero00.grid(column=0, row=20, sticky = 'WES')
    vari = ['True', 'False']
    mostrar_vari = StringVar()
    vari_in_chord = ttk.Combobox(root, textvariable = mostrar_vari, font=("Arial", SIZE_TEXT_ALL),
                         values = vari, width=ancho_boton)
    vari_in_chord.grid(column=1, row=20, sticky= 'SW')
    vari_in_chord.bind('<<ComboboxSelected>>', on_select)
    vari_in_chord.current(0)


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    cero1 = Label(root, text='Save plot', font=("Arial", SIZE_TEXT_ALL, 'bold'), fg = 'black', bg = 'lime')
    cero1.grid(column=0, row=21, sticky = 'WES')




    #-------------------------------------
    #-------------------------------------


    agrupar = LabelFrame(root, text = "", font=("Arial", 8))
    agrupar.grid(column=0, row=22, columnspan = 2, sticky= 'WENS')
    agrupar.configure(background='white')


    def on_select2():
        fig.savefig('plots_asv/taxonomy/Chord_'+title_kit+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.png', dpi = 900, bbox_inches= 'tight')

    ###############
    ###############


    boton = Button(agrupar, text="PNG", cursor="hand2", borderwidth=0,
                bg="gold", fg="black",font=("Arial", SIZE_TEXT_ALL), command = on_select2)
    boton.grid(column = 0, row = 0, sticky = 'WES')
    #boton.bind('<Button-1>', on_select2)
    
    labebl = Label(agrupar, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 1, row = 0)

    def on_select3():
        fig.savefig('plots_asv/taxonomy/Chord_'+title_kit+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.jpeg', dpi = 900, bbox_inches= 'tight')
    ###############
    ###############


    boton = Button(agrupar, text="JPEG", cursor="hand2", borderwidth=0,
                bg="gold", fg="black",font=("Arial", SIZE_TEXT_ALL), command = on_select3)
    boton.grid(column = 2, row = 0, sticky = 'WES')
    #boton.bind('<Button-1>', on_select2)
    
    labebl = Label(agrupar, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 3, row = 0)

    #-------------------------------------
    #-------------------------------------
    def on_select4():
        fig.savefig('plots_asv/taxonomy/Chord_'+title_kit+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.svg', dpi = 900, bbox_inches= 'tight')
    ###############
    ###############


    boton = Button(agrupar, text="SVG", cursor="hand2", borderwidth=0,
                bg="gold", fg="black",font=("Arial", SIZE_TEXT_ALL), command = on_select4)
    boton.grid(column = 4, row = 0, sticky = 'WES')
    #boton.bind('<Button-1>', on_select3)
    
    labebl = Label(agrupar, text= '', font=("Arial", 8), fg="red", bg = 'white')
    labebl.grid(column = 5, row = 0)
    #-------------------------------------
    #-------------------------------------
    def on_select5():
        fig.savefig('plots_asv/taxonomy/Chord_'+title_kit+'_'+Linaje+'_'+Porcentaje+'_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.pdf', dpi = 900, bbox_inches= 'tight')
    ###############
    ###############


    boton = Button(agrupar, text="PDF", cursor="hand2", borderwidth=0,
                bg="gold", fg="black",font=("Arial", SIZE_TEXT_ALL), command = on_select5)
    boton.grid(column = 6, row = 0, sticky = 'WES')
    #boton.bind('<Button-1>', on_select3)


    #------------------------------------


    group_aspect = LabelFrame(root, text = ' CHORT PLOT:  '+title_kit+' ', font=("Arial", 10,  "bold"), height = 5, fg = 'black', bg = 'blue')
    group_aspect.grid(column=2, row=0, rowspan = 200, sticky= 'EN')
    group_aspect.configure(background='white')

    mpl.rcParams.update(mpl.rcParamsDefault)

    fig = plt.figure(figsize=(7, 7))

    bar1 = FigureCanvasTkAgg(fig, group_aspect)
    bar1.draw()
    bar1.get_tk_widget().grid(row=0, column = 0)#, rowspan=7, columnspan = 7


    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_aspect('equal', 'box')
    #ax.set_facecolor('none')
    

    sources_label = 1.36
    targets_label = 1.58

    radio = float(radio_fuente[9])
    sepp = 0.01
    ancho = ancho_sources[7] # % del ancho
    #tam = 5
    Espacio = 0.5 # interespacio_source[5]


    constante = apertura[160] - len(origenes)
    tam = constante / len(origenes) # 

    W = tam * len(origenes)

    Q = (Espacio * len(origenes)) - Espacio
    sour_inicial = 360 - ((W + Q)/ 2)
    #sour_inicial = 360 - 90

    """
    Nodos de Sources sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
    """
    nodos_ori_0 = {}
    dif0 = 360 - ((W + Q)/ 2)
    dif1 = dif0 # inicio de los nodos
    record_ori_0 = {}

    for o in origenes:
        teta1 = dif0
        teta2 = dif0+tam

        sour, mitad3, etiquetapos, mitad1 = source1(radio = radio, theta1 = teta1, theta2 = teta2,
                                                 width = radio*(ancho/100), sep = sepp, color = CorrepondenciA[o])
        ax.add_patch(sour)
        record_ori_0[o] = [teta1, teta2]

        central_angle0 = (((dif0+tam)-dif0)/2)+dif0
        central_angle = (central_angle0*np.pi)/180
        if etiquetapos[0] < 0:
            central_angle = central_angle - np.pi
            if len(o) < pal_max_origenes:
                tam_pal = len(o)
                palabra = ' '*(pal_max_origenes - tam_pal)+o
                palabra = palabra+' '*len(palabra)+' '
            else:
                palabra = o
                palabra = palabra+' '*len(palabra)+' '
        else:
            central_angle = (central_angle0*np.pi)/180
            if len(o) < pal_max_origenes:
                tam_pal = len(o)
                palabra = o+' '*(pal_max_origenes - tam_pal)
                palabra = ' '*len(palabra)+' '+palabra
            else:
                palabra = o
                palabra = ' '*len(palabra)+' '+palabra


        ax.text(etiquetapos[0], etiquetapos[1], palabra, color = 'black',
                va = 'center', ha = 'center', #fontweight = 'bold',
                fontsize = letra_source[6], rotation = np.rad2deg(central_angle),
                family = 'monospace', style='italic')
        #ax.scatter(etiquetapos[0], etiquetapos[1], s = 50, c = 'black', zorder = 2)

        if Counts_SOURCE[o] == 1:
            t1 = dif1
            t2 = dif1+tam
            sour2, mitad33, internopos = source2(radio = radio, theta1 = t1, theta2 = t2,
                                                 width = ancho, sep = sepp, color = 'white')
            ax.add_patch(sour2)
            nodos_ori_0[o] = [sour2, mitad33, internopos]
            dif1 += (tam+Espacio)
        else:
            t1 = dif1
            sectores = tam/Counts_SOURCE[o]
            ss0 = 0
            ss1 = sectores
            intersecciones_ori = []
            for r in range(Counts_SOURCE[o]):
                t1 = ss0+dif1
                t2 = ss1+dif1
                sour2, mitad33, internopos = source2(radio = radio, theta1 = t1, theta2 = t2,
                                                     width = ancho, sep = sepp, color = 'white')

                ax.add_patch(sour2)
                intersecciones_ori.append([sour2, mitad33, internopos])

                ss1 += sectores
                ss0 += sectores
            nodos_ori_0[o] = intersecciones_ori
            dif1 += (tam+Espacio)

        PER = 2
        central_angle = (((dif0+tam)-dif0)/2)+dif0

        dif0 += (tam+Espacio)

    """
    Nodos de Targets ttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
    """   
    #----------
    # angulo de separacion elegido

    angulo_espacio = tamanos[45] # tamanos[55]

    continuacion = dif0 - Espacio - 360
    Espacio2 = 0.5

    if angulo_espacio > 0:
        dif00 = angulo_espacio + continuacion
        tam2 = ((sour_inicial-angulo_espacio)-(continuacion+angulo_espacio) - ((len(destinos) * Espacio2)-Espacio2))/len(destinos)
    else:
        tam2 = ((sour_inicial-continuacion) - ((len(destinos) * Espacio2)-Espacio2))/len(destinos)
        dif00 = continuacion

    nodos_des_0 = {}

    ANGULOS = []
    nodos_des_1 = {}
    #dif00 = continuacion
    dif11 = dif00
    #tam2 = ((sour_inicial-continuacion) - ((len(destinos) * Espacio2)-Espacio2))/len(destinos)

    radio2 = float(radio_destino[9])

    sepp2 = 0.01
    ancho2 = ancho_targets[5]

    var_radios = {list(dict_variable_element_colors.keys())[0]:radio2}



    increase = sepp2+(radio2*(ancho2/100))
    for i in list(dict_variable_element_colors.keys())[1:]:
        var_radios.update({i:radio2+increase})
        increase += sepp2+(radio2*(ancho2/100))

    MitadeS = []

    for o in destinos:
        teta11 = dif11
        teta22 = dif11+tam2

        if 'None' in list(dict_variable_element_colors.keys()):

            sour, mitad2, etiquetapos, mitad1 = source1(radio = radio2, theta1 = teta11, theta2 = teta22,
                                                     width = radio2*(ancho2/100), sep = sepp2, color = 'silver')
            ax.add_patch(sour)
            ETIQUETAS = etiquetapos

            MitadeS.append([mitad1, mitad2])

        else:
            #### Location
            if 'Location' in list(dict_variable_element_colors.keys()):
                sour, mitad2, etiquetapos, mitad1 = source1(radio = var_radios['Location'], theta1 = teta11, theta2 = teta22,
                                                         width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Location'][correspondencia_sam_vars[o]['Location']])
                ax.add_patch(sour)
                ETIQUETAS = etiquetapos

                if o == destinos[0]:
                    MITADES = [mitad1, mitad2]
                    MitadeS.append(MITADES)

            ##### ug OTA/kg
            if 'ug OTA/kg' in list(dict_variable_element_colors.keys()):
                sour0, mitad22, etiquetapos0, mitad11 = source1(radio = var_radios['ug OTA/kg'], theta1 = teta11, theta2 = teta22,
                                                         width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['ug OTA/kg'][correspondencia_sam_vars[o]['ug OTA/kg']])
                ax.add_patch(sour0)
                ETIQUETAS = etiquetapos0

                if o == destinos[0]:
                    MITADES = [mitad11, mitad22]
                    MitadeS.append(MITADES)


            ##### Cultivation
            if 'Cultivation' in list(dict_variable_element_colors.keys()):
                sour00, mitad222, etiquetapos00, mitad111 = source1(radio = var_radios['Cultivation'], theta1 = teta11, theta2 = teta22,
                                                         width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Cultivation'][correspondencia_sam_vars[o]['Cultivation']])
                ax.add_patch(sour00)
                ETIQUETAS = etiquetapos00

                if o == destinos[0]:
                    MITADES = [mitad111, mitad222]
                    MitadeS.append(MITADES)

            ##### Coffee Variety
            if 'Coffee Variety' in list(dict_variable_element_colors.keys()):
                sour000, mitad2222, etiquetapos000, mitad1111 = source1(radio = var_radios['Coffee Variety'], theta1 = teta11, theta2 = teta22,
                                                         width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Coffee Variety'][correspondencia_sam_vars[o]['Coffee Variety']])
                ax.add_patch(sour000)
                ETIQUETAS = etiquetapos000

                if o == destinos[0]:
                    MITADES = [mitad1111, mitad2222]
                    MitadeS.append(MITADES)

            ##### Genomic DNA kit
            if 'Genomic DNA kit' in list(dict_variable_element_colors.keys()):
                sour0000, mitad22222, etiquetapos0000, mitad11111 = source1(radio = var_radios['Genomic DNA kit'], theta1 = teta11, theta2 = teta22,
                                                         width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Genomic DNA kit'][correspondencia_sam_vars[o]['Genomic DNA kit']])
                ax.add_patch(sour0000)
                ETIQUETAS = etiquetapos0000

                if o == destinos[0]:
                    MITADES = [mitad11111, mitad22222]
                    MitadeS.append(MITADES)


            ##### Drying Time (Days)
            if 'Drying Time (Days)' in list(dict_variable_element_colors.keys()):
                sour00000, mitad222222, etiquetapos00000, mitad111111 = source1(radio = var_radios['Drying Time (Days)'], theta1 = teta11, theta2 = teta22,
                                                         width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Drying Time (Days)'][correspondencia_sam_vars[o]['Drying Time (Days)']])
                ax.add_patch(sour00000)
                ETIQUETAS = etiquetapos00000

                if o == destinos[0]:
                    MITADES = [mitad111111, mitad222222]
                    MitadeS.append(MITADES)

            ##### Postharvest Processing
            if 'Postharvest Processing' in list(dict_variable_element_colors.keys()):
                sour000000, mitad2222222, etiquetapos000000, mitad1111111 = source1(radio = var_radios['Postharvest Processing'], theta1 = teta11, theta2 = teta22,
                                                         width = radio2*(ancho2/100), sep = sepp2, color = dict_variable_element_colors['Postharvest Processing'][correspondencia_sam_vars[o]['Postharvest Processing']])
                ax.add_patch(sour000000)
                ETIQUETAS = etiquetapos000000

                if o == destinos[0]:
                    MITADES = [mitad1111111, mitad2222222]
                    MitadeS.append(MITADES)

        #####
        central_angle0 = (((dif00+tam2)-dif00)/2)+dif00
        central_angle = (central_angle0*np.pi)/180


        if etiquetapos[0] < 0:
            central_angle = central_angle - np.pi
            if len(o) < pal_max_destinos:
                tam_pal = len(o)
                palabra = ' '*(pal_max_destinos - tam_pal)+o
                palabra = palabra+' '+' '*len(palabra)
            else:
                palabra = o
                palabra = palabra+' '+' '*len(palabra)

        else:
            central_angle = (central_angle0*np.pi)/180
            if len(o) < pal_max_destinos:
                tam_pal = len(o)
                palabra = o+' '*(pal_max_destinos - tam_pal)
                palabra = ' '+' '*len(palabra)+palabra
            else:
                palabra = o
                palabra = ' '+' '*len(palabra)+palabra


        ANGULOS.append(dif00)
        ax.text(ETIQUETAS[0], ETIQUETAS[1], palabra, color = 'black', va = 'center', ha = 'center', 
                fontsize = letra_target[5], rotation = np.rad2deg(central_angle), fontfamily = 'Liberation mono', fontweight='bold')

        #ax.scatter(etiquetapos33[0], etiquetapos33[1], s = 20, c = 'lime', zorder = 2)

        if Counts_TARGET[o] == 1:
            t1 = dif11
            t2 = dif11+tam2
            sour2, mitad33, internopos = source2(radio = radio2, theta1 = t1, theta2 = t2,
                                                 width = ancho2, sep = sepp2, color = 'white')
            ax.add_patch(sour2)
            nodos_des_1[o] = [sour2, mitad33, internopos]
            dif11 += (tam2+Espacio2)
        else:
            t1 = dif11
            sectores = tam2/Counts_TARGET[o]
            ss0 = 0
            ss1 = sectores
            intersecciones_des = []
            for r in range(Counts_TARGET[o]):
                t1 = ss0+dif11
                t2 = ss1+dif11
                sour2, mitad33, internopos = source2(radio = radio2, theta1 = t1, theta2 = t2,
                                                     width = ancho2, sep = sepp2, color = 'white')
                ax.add_patch(sour2)
                intersecciones_des.append([sour2, mitad33, internopos])
                ss1 += sectores
                ss0 += sectores
            nodos_des_1[o] = intersecciones_des
            dif11 += (tam2+Espacio2)

        central_angle2 = (((dif00+tam2)-dif00)/2)+dif00
        radian2 = (central_angle2*np.pi)/180
        x2 = (radio2 * (1-((ancho2/2)/100))) *  np.cos(radian2)
        y2 = (radio2 * (1-((ancho2/2)/100))) *  np.sin(radian2)
        #ax.scatter(x2, y2, s = 10, c = 'black', zorder = 2)

        dif00 += (tam2+Espacio2)

    """
    Conexiones ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    """
    #---------------
    # parte 1
    XX = []
    for ori in origenes:
        xu = 0
        if Counts_SOURCE[ori] == 1:
            for SOURCE, TARGET in net:
                if ori == SOURCE:
                    if Counts_TARGET[TARGET] == 1:  # target uno
                        #print('>>>', SOURCE, TARGET)
                        path_data = LOCATIONS(nodos_ori_0[SOURCE][1], nodos_des_1[TARGET][1])

                        codes, verts = zip(*path_data)
                        path = mpath.Path(verts, codes)
                        patch = mpatches.PathPatch(path, facecolor=CorrepondenciA[ori], alpha=1, lw = None, ec = 'none', zorder = 0)
                        ax.add_patch(patch)
                    else: # target con mas de uno
                        XX.append([SOURCE, 'NA', TARGET])
        else:
            for SOURCE, TARGET in net:
                if ori == SOURCE:
                    #print(SOURCE, xu, TARGET)
                    if Counts_TARGET[TARGET] == 1:
                        #print(SOURCE, xu, TARGET, '----')
                        path_data = LOCATIONS(nodos_ori_0[SOURCE][xu][1], nodos_des_1[TARGET][1])
                        codes, verts = zip(*path_data)
                        path = mpath.Path(verts, codes)
                        patch = mpatches.PathPatch(path, facecolor=CorrepondenciA[ori], alpha=1, lw = None, ec = 'none', zorder = 0)
                        ax.add_patch(patch)
                    else: # target con mas de uno
                        XX.append([SOURCE, xu, TARGET])
                    xu += 1


    #---------------
    # Parte 2
    output = []
    for SOURCE, P, TARGET in XX:
        if TARGET not in output:
            output.append(TARGET)

    for s in output:
        n = 0
        for SOURCE, P, TARGET in XX:
            if s == TARGET:

                if P == 'NA':
                    path_data = LOCATIONS(nodos_ori_0[SOURCE][1], nodos_des_1[TARGET][n][1])

                    codes, verts = zip(*path_data)
                    path = mpath.Path(verts, codes)
                    patch = mpatches.PathPatch(path, facecolor=CorrepondenciA[SOURCE], alpha=1, lw = None, ec = 'none', zorder = 0)
                    ax.add_patch(patch)
                else:
                    path_data = LOCATIONS(nodos_ori_0[SOURCE][P][1], nodos_des_1[TARGET][n][1])

                    codes, verts = zip(*path_data)
                    path = mpath.Path(verts, codes)
                    patch = mpatches.PathPatch(path, facecolor=CorrepondenciA[SOURCE], alpha=1, lw = None, ec = 'none', zorder = 0)
                    ax.add_patch(patch)
                n += 1
    #

    for i, j in zip(MitadeS, list(dict_variable_element_colors.keys())):
        ax.text(np.sum(i[0][0][0]+i[1][0][0])/2, np.sum(i[0][0][1]+i[1][0][1])/2,
                ' '*len(j)+' '+j, color = 'black', ha = 'center',va = 'center', #weight = 'bold',
                fontsize = 6, rotation = ((continuacion+angulo_espacio)-90)-2, fontfamily = 'Liberation mono')

    #
    #ax.text(0, 1.9, title_kit, ha = 'center', va = 'center', fontsize = 10, weight = 'bold')


    if 'None' in list(dict_variable_element_colors.keys()):
        pass
    else:
        yy = 1.8
        xx = -2.4
        for var in list(dict_variable_element_colors.keys()):
            ax.text(xx-0.09, yy, var, ha = 'left', va = 'center', fontsize = 7, fontfamily = 'Arial')
            yy -= 0.07
            for e, j in enumerate(dict_variable_element_colors[var]):
                ax.scatter(xx-0.07, yy, s = 10, c = dict_variable_element_colors[var][j], marker = 'o')
                ax.text(xx-0.05, yy, '  '+j, ha = 'left', va = 'center', fontsize = 7, rotation = 0, fontfamily = '3ds Light')
                #xx +=0.07
                yy -= 0.07
            yy -= 0.05


    ax.set_xlim(-float(limites[1])-0.5, float(limites[1]))
    ax.set_ylim(-float(limites[1]), float(limites[1]))
    ax.axis('off')

    plt.close()


    root.mainloop()





Chord = widgets.Button(description="CHORT PLOT", icon = 'fa-circle-o', layout=Layout(width='120px', height='27px'))
Chord.style.button_color = 'gold'
Chord_out = widgets.Output()
def button_clicked1(b):
    with Chord_out:
        clear_output(True)
        chord_plot()
Chord.on_click(button_clicked1)











metricas = ['euclidean','braycurtis','canberra','chebyshev','cityblock','correlation','cosine','dice',
            'hamming','jaccard','jensenshannon','kulsinski','mahalanobis',
            'matching','minkowski','rogerstanimoto','russellrao','seuclidean','sokalmichener',
            'sokalsneath','sqeuclidean','yule']
metodos = ['complete','single','average','weighted','centroid','ward']








uno = widgets.HTML('<i class="fa fa-spinner fa-2x fa-fw"></i>')


TAXONOMY_button = widgets.Button(description="PROCESS AND VISUALIZE", icon = 'fa-eye', layout=Layout(width='590px'))
TAXONOMY_button.style.button_color = 'gainsboro' #'deepskyblue'
TAXONOMY_button.style.font_weight = 'bold'
TAXONOMY_output = widgets.Output()

va = []
n = 0.1
for i in range(150):
    va.append(round(n, 3))
    n += 0.1
va = [0, 0.01, 0.05] + va


def button_clicked(b):
    
    from scipy.cluster.hierarchy import dendrogram, linkage
    from matplotlib.colors import ListedColormap # funcion para crear un objeto <matplotlib.colors.ListedColormap> a partir de una lista de colores personalizados
    
    
    NCBI_RDP_SILVA_SUMMARY = pd.read_csv('tablas/ASVs_NCBI_RDP_SILVA_SUMMARY.txt', sep = '\t')

    ASV_Full_Taxonomy = pd.read_csv('tablas/ASV_Full_Taxonomy.txt', sep = '\t')
    ASV_Full_Taxonomy = ASV_Full_Taxonomy.merge(ASV_counts, on = 'Entry', how = 'left')
    
    with TAXONOMY_output:
        clear_output(True)
        
        tipo_kits_1 = widgets.ToggleButtons(options= ['Both kits'] + metadata[VARIABLE_KIT].unique().tolist(), value = 'Both kits', button_style = 'primary')
        tipo_kits_1.style.button_width = '170px'
        tipo_kits_1.style.font_weight = 'bold'
        tipo_kit_1_box = Box(children=[VBox([tipo_kits_1])], layout=Layout(border='1px solid #1976d2', width='180px', height='95px'))
        
        tax = widgets.ToggleButtons(options= category_names, value = 'Genus', button_style = '')
        tax.style.button_width = '80px'
        tipo_tax_box = Box(children=[VBox([tax])], layout=Layout(border='1px solid gainsboro', width='90px', height='210px'))
        
    
        PercentagE = widgets.SelectionSlider(options=list(reversed(va)),value=1,disabled=False, continuous_update=False,
                             orientation='vertical',readout=True, layout=Layout(width='25px', height='180px'))
        PercentagE_box = VBox([widgets.HTML('<b>Limit (%):</b>'), PercentagE])
        
        METODO = widgets.Dropdown(options=metodos,value='complete',disabled=False,
                                 layout = Layout(width='290px', height='25px'))
        METRICA = widgets.Dropdown(options=metricas,value='euclidean',disabled=False,
                                 layout = Layout(width='290px', height='25px'))
        
        para_mostrar = widgets.ToggleButtons(options=['Top 10', 'Top 15', 'Top 20', 'All'], value = 'All')
        para_mostrar.style.button_width = '70px'
        
        
        Multiple_colorS = widgets.SelectMultiple(options=sorted(list(qualitative_colors.keys())), value=sorted(list(qualitative_colors.keys()))[0:7], disabled=False,
                       layout=Layout(width='120px', height='250px'))
        
        
        tam_plot1 = widgets.SelectionSlider(options=[9, 10, 11, 12],value=10,disabled=False,
                                              description = 'Chart size:', continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        n = 0
        wl = []
        for e, i in enumerate(range(21)):
            wl.append(round(n, 2))
            n += 0.2
        width_linea = widgets.SelectionSlider(options=wl, value=1.4,disabled=False,
                                              description = 'Line width:', continuous_update=False,orientation='horizontal',readout=True,
                                             layout=Layout(width='400px', height='25px'))
        
        family = sorted(['Liberation Serif','Microsoft Sans Serif','Open Sans','Times New Roman','3ds Light','Calibri','Comic Sans MS',
                  'Arial','Courier New','Microsoft Yi Baiti','Lucida Console'])

        
        family_axis = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                   layout = Layout(width='290px', height='25px'))
        
        despazamiento = {}
        n = 1.1
        for e, i in enumerate(range(70)):
            despazamiento[e+1] = round(n, 3)
            n += 0.01

        displacement_leye = widgets.SelectionSlider(options=list(despazamiento.keys()),value= 36,disabled=False,
                                                      description = 'Legend pos:',
                                                continuous_update=False,orientation='horizontal',readout=True,
                                                   layout=Layout(width='400px', height='25px'))
        
        num_cols = widgets.ToggleButtons(options=[1, 2, 3, 4, 5, 5, 7, 8], value = 1)
        num_cols.style.button_width = '30px'
        
        
        
        Variables_metadatA = widgets.SelectMultiple(options=['None'] + variables, value=variables, disabled=False, layout = Layout(width='180px', height='180px'))
        Variables_metadatA_box = Box(children=[VBox([Variables_metadatA])], layout=Layout(border='2px solid #097cd7', width='187px', height='187px'))
        
        sample_text = widgets.SelectionSlider(options=np.arange(5, 12.5, 0.5).tolist(),value=7.5,disabled=False,
                                              description = 'Sample text:', continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        
        varia_text = widgets.SelectionSlider(options=np.arange(5, 12.5, 0.5).tolist(),value=6.5,disabled=False,
                                              description = 'Variable text:', continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        
        ancho_dendo = widgets.SelectionSlider(options=[0.05, 0.075, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2],value= 0.15,disabled=False,
                                              description = 'Dendrogram width:', continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
         
        cambiar_sample = widgets.ToggleButtons(options=['Code1', 'Code2', 'Code3'], value = 'Code1')
        cambiar_sample.style.button_width = '70px'
        
        ######

        VAR0 = 0
        VAR1 = 1
        VAR2 = 2
        VAR3 = 3
        VAR4 = 4
        VAR5 = 5
        VAR6 = 6

        width_cuadro = str(100)

        blanca3 = widgets.Button(layout=Layout(width='10px', height='2px'), disabled=True)
        blanca3.style.button_color = 'white'

        for unoo in [variables[VAR0]]:
            Lab_Var0 = widgets.Label(unoo)
            lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
            SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
            def show_rampa(SalidA0):
                #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                patron = SalidA0.split('_v')[0]
                barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
            OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})

        for dos in [variables[VAR1]]:
            Lab_Var1 = widgets.Label(dos)
            lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[dos])].keys())
            SalidA1 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
            def show_rampa(SalidA1):
                #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[dos])][SalidA1])
                patron = SalidA1.split('_v')[0]
                barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[dos])][SalidA1])
            OUTshow_rampa1 = widgets.interactive_output(show_rampa, {'SalidA1':SalidA1})

        for tres in [variables[VAR2]]:
            Lab_Var2 = widgets.Label(tres)
            lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[tres])].keys())
            SalidA2 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
            def show_rampa(SalidA2):
                #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[tres])][SalidA2])
                patron = SalidA2.split('_v')[0]
                barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[tres])][SalidA2])
            OUTshow_rampa2 = widgets.interactive_output(show_rampa, {'SalidA2':SalidA2})

        for cuatro in [variables[VAR3]]:
            Lab_Var3 = widgets.Label(cuatro)
            lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[cuatro])].keys())
            SalidA3 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
            def show_rampa(SalidA3):
                #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[cuatro])][SalidA3])
                patron = SalidA3.split('_v')[0]
                barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[cuatro])][SalidA3])
            OUTshow_rampa3 = widgets.interactive_output(show_rampa, {'SalidA3':SalidA3})

        for cinco in [variables[VAR4]]:
            Lab_Var4 = widgets.Label(cinco)
            lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[cinco])].keys())
            SalidA4 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
            def show_rampa(SalidA4):
                #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[cinco])][SalidA4])
                patron = SalidA4.split('_v')[0]
                barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[cinco])][SalidA4])
            OUTshow_rampa4 = widgets.interactive_output(show_rampa, {'SalidA4':SalidA4})

        for seis in [variables[VAR5]]:
            Lab_Var5 = widgets.Label(seis[0:13] + ' ...')
            lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[seis])].keys())
            SalidA5 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
            def show_rampa(SalidA5):
                #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[seis])][SalidA5])
                patron = SalidA5.split('_v')[0]
                barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[seis])][SalidA5])
            OUTshow_rampa5 = widgets.interactive_output(show_rampa, {'SalidA5':SalidA5})

        for siete in [variables[VAR6]]:
            Lab_Var6 = widgets.Label(siete[0:13] + ' ...')
            lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[siete])].keys())
            SalidA6 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
            def show_rampa(SalidA6):
                #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[siete])][SalidA6])
                patron = SalidA6.split('_v')[0]
                barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[siete])][SalidA6])
            OUTshow_rampa6 = widgets.interactive_output(show_rampa, {'SalidA6':SalidA6})

        Wid_Vars = [HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0]),
                    HBox([VBox([blanca3, VBox([Lab_Var1, SalidA1])]), OUTshow_rampa1]),
                    HBox([VBox([blanca3, VBox([Lab_Var2, SalidA2])]), OUTshow_rampa2]),
                    HBox([VBox([blanca3, VBox([Lab_Var3, SalidA3])]), OUTshow_rampa3]),
                    HBox([VBox([blanca3, VBox([Lab_Var4, SalidA4])]), OUTshow_rampa4]),
                    HBox([VBox([blanca3, VBox([Lab_Var5, SalidA5])]), OUTshow_rampa5]),
                    HBox([VBox([blanca3, VBox([Lab_Var6, SalidA6])]), OUTshow_rampa6])]

        Colores_para_variableS = Box(children = Wid_Vars, layout=Layout(flex_flow='column', align_items='flex-start', border='1px solid gainsboro', width='280px', height='540px')) # 
        
        filename = widgets.Text(value='', placeholder='File name.txt', description='', disabled=False, layout = Layout(width='195px', height='25px'))
        
        filename_plot = widgets.Text(value='', placeholder='Chart name', description='', disabled=False, layout = Layout(width='270px', height='25px'))
        ######
        
        
        def TaxonomY(tipo_kits_1, tax, PercentagE, METODO, METRICA, para_mostrar, Multiple_colorS, tam_plot1, width_linea, Variables_metadatA, family_axis, displacement_leye,
                    num_cols, SalidA0, SalidA1, SalidA2, SalidA3, SalidA4, SalidA5, SalidA6, sample_text, varia_text, ancho_dendo, cambiar_sample):
            
            #uno.value = '<i class="fa fa-spinner fa-pulse fa-2x fa-fw"></i> <font color = black>  <h style="font-size:0.6vw"> Processing.</font>'
            uno.value = '<i class="fa fa-spinner fa-pulse fa-2x fa-fw"></i>'

            

            if tipo_kits_1 == 'Both kits':
                rarefaction_collapse_specie = pd.pivot_table(NCBI_RDP_SILVA_SUMMARY[['Species'] + name_sample], values = name_sample, index = ['Species'], aggfunc = sum).reset_index()
                rarefaction_collapse_specie['Sum'] = np.sum(rarefaction_collapse_specie[name_sample].values, axis = 1)
                NEW_COLUMNS = name_sample

            if tipo_kits_1 == 'DNPowerSoil':
                rarefaction_collapse_specie = pd.pivot_table(NCBI_RDP_SILVA_SUMMARY[['Species'] + DNPowerSoil], values = DNPowerSoil, index = ['Species'], aggfunc = sum).reset_index()
                rarefaction_collapse_specie['Sum'] = np.sum(rarefaction_collapse_specie[DNPowerSoil].values, axis = 1)
                rarefaction_collapse_specie = rarefaction_collapse_specie[rarefaction_collapse_specie.Sum > 0]
                rarefaction_collapse_specie = rarefaction_collapse_specie.drop(columns = ['Sum'])
                NEW_COLUMNS = DNPowerSoil

            if tipo_kits_1 == 'DNMicrobial':
                rarefaction_collapse_specie = pd.pivot_table(NCBI_RDP_SILVA_SUMMARY[['Species'] + DNMicrobial], values = DNMicrobial, index = ['Species'], aggfunc = sum).reset_index()
                rarefaction_collapse_specie['Sum'] = np.sum(rarefaction_collapse_specie[DNMicrobial].values, axis = 1)
                rarefaction_collapse_specie = rarefaction_collapse_specie[rarefaction_collapse_specie.Sum > 0]
                rarefaction_collapse_specie = rarefaction_collapse_specie.drop(columns = ['Sum'])
                NEW_COLUMNS = DNMicrobial
                
            counts_lineaje = rarefaction_collapse_specie.merge(ASV_Full_Taxonomy[category_names].drop_duplicates(), on = 'Species', how = 'left')
            data_pivot = pd.pivot_table(counts_lineaje[[tax] + NEW_COLUMNS], values = NEW_COLUMNS, index = [tax], aggfunc = sum).reset_index()

            percentages = data_pivot[NEW_COLUMNS].values/data_pivot[NEW_COLUMNS].values.sum(axis=0)*100 # porcentajes
            relative_abundance = DataFrame(percentages, columns = NEW_COLUMNS)
            relative_abundance.insert(loc = 0, column=tax, value=data_pivot[tax])

            abundancias = []
            for i in NEW_COLUMNS:
                df = relative_abundance[[tax, i]]
                if PercentagE == 0:
                    df = df[df[i] > PercentagE]
                    ff = df
                else:
                    df = df[df[i] >= PercentagE]
                    ef = DataFrame({tax:['Others'], i:[100 - df[i].sum()]})
                    ff = pd.concat([df, ef])
                abundancias.append(ff)
            AbundancE = reduce(lambda  left,right: pd.merge(left, right, on = [tax], how = 'outer'), abundancias).fillna(0)

            

            # metodo y metrica
            #METODO = 'complete'
            #METRICA = 'euclidean'
            
            dendro = dendrogram(linkage(AbundancE[NEW_COLUMNS].values.T, method = METODO, metric = METRICA),
                        orientation='left',  no_plot = True,
                        labels=AbundancE[NEW_COLUMNS].columns.tolist(),
                        distance_sort='descending',
                                leaf_font_size = 100,
                        show_leaf_counts=True)
            #-----------------------------------------
            ivl = dendro['ivl']
            color_list = dendro['color_list']
            Z = np.asarray(dendro['dcoord'], order='c')
            mh = max(Z[:, 2])
            above_threshold_color = 'b'
            # Independent variable plot width
            ivw = len(ivl) * 10
            # Dependent variable plot height
            dvw = mh + mh * 0.05
            iv_ticks = np.arange(5, len(ivl) * 10 + 5, 10)
            xlines = dendro['dcoord']
            ylines = dendro['icoord']
            
            
            AAAbundancEEE = AbundancE[dendro['ivl']].T
            AAAbundancEEE.columns = AbundancE[tax].tolist()

            # ordena las columnas de mayor a menor
            AAAbundancEEE = AAAbundancEEE[list(dict(OrderedDict(Counter(dict(zip(AAAbundancEEE.columns, np.sum(AAAbundancEEE.values, axis = 0)))).most_common())).keys())]
            
            
            
            ColS = AAAbundancEEE.columns.tolist()
            
            if para_mostrar == 'All':
                AAAbundancEEE.to_csv('plots_asv/taxonomy/AAAbundancEEE.txt', sep = '\t')
                pass
            if para_mostrar == 'Top 10':
                if PercentagE == 0:
                    AAAbundancEEE = AAAbundancEEE[ColS[:10]]
                    AAAbundancEEE.to_csv('plots_asv/taxonomy/AAAbundancEEE.txt', sep = '\t')
                else:
                    ColS.remove('Others')
                    AAAbundancEEE = AAAbundancEEE[ColS[:10]]
                    AAAbundancEEE.to_csv('plots_asv/taxonomy/AAAbundancEEE.txt', sep = '\t')
            if para_mostrar == 'Top 15':
                if PercentagE == 0:
                    AAAbundancEEE = AAAbundancEEE[ColS[:15]]
                    AAAbundancEEE.to_csv('plots_asv/taxonomy/AAAbundancEEE.txt', sep = '\t')
                else:
                    ColS.remove('Others')
                    AAAbundancEEE = AAAbundancEEE[ColS[:15]]
                    AAAbundancEEE.to_csv('plots_asv/taxonomy/AAAbundancEEE.txt', sep = '\t')
            if para_mostrar == 'Top 20':
                if PercentagE == 0:
                    AAAbundancEEE = AAAbundancEEE[ColS[:20]]
                    AAAbundancEEE.to_csv('plots_asv/taxonomy/AAAbundancEEE.txt', sep = '\t')
                else:
                    ColS.remove('Others')
                    AAAbundancEEE = AAAbundancEEE[ColS[:20]]
                    AAAbundancEEE.to_csv('plots_asv/taxonomy/AAAbundancEEE.txt', sep = '\t')
                
            rangonumeros = np.arange(len(AAAbundancEEE.columns.tolist()))
            serie = dict(zip(AAAbundancEEE.columns.tolist(), rangonumeros))


            # rampa = 'tab20c' o ListedColormap(Set1+Set2+Pastel1+Pastel2+tab20+tab20c+Paired)
            """
            creacion de rampa de colores sin duplicados
            """
            rampa_creada = list(dict.fromkeys([to_hex(k) for i in Multiple_colorS for k in plt.get_cmap(i)(np.arange(qualitative_colors[i])/qualitative_colors[i])]))
            
            cuenta_colors_tax = HBox([widgets.Label('BAR COLORS='+str(len(rampa_creada))+', '+tax+'='+str(len(serie)))])
            
            if len(serie) <= len(rampa_creada):
                rampa_creada = rampa_creada[0:len(serie)]
            else:
                pass
            
            
            
            cnorm = mpl.colors.Normalize(vmin=0, vmax=max(list(serie.values())))
            cpick = cm.ScalarMappable(norm=cnorm, cmap= ListedColormap(rampa_creada))

            cpick.set_array([])
            val_map = {}
            for k, v in zip(list(serie.keys()), list(serie.values())):
                val_map[k] = cpick.to_rgba(v)
            colors = [] # rgb
            colors2 = {} # hex
            dict_para_salvar = {}
            for node in list(serie.keys()):
                colors.append(val_map[node])
                colors2[node] = to_hex(val_map[node])
                dict_para_salvar[node] = to_hex(val_map[node])
            
            dict_para_salvar.update({'data':tipo_kits_1, 'Linaje':tax, 'Porcentaje':str(PercentagE)})
            with open('plots_asv/taxonomy/dict_para_salvar.json', 'w') as fp:
                json.dump(dict_para_salvar, fp)
                
            
            mues = {}
            for indi in AAAbundancEEE.index:
                df = AAAbundancEEE[AAAbundancEEE.index == indi]
                df= df[df > 0].dropna(axis=1)
                """
                ordenada
                """
                mues[indi] = dict(OrderedDict(Counter(dict(zip(df.columns, df.values[0]))).most_common()))
                """
                sin modificar
                """
                #mues[indi] = dict(zip(df.columns, df.values[0]))
                
            ### asignacion de colores a cada variable, se asigna una lista
            
            asignaciones_colors_variables = {'Location':diccionarios_colores[len(variables_items_ordened['Location'])][SalidA0],
                                             'ug OTA/kg':diccionarios_colores[len(variables_items_ordened['ug OTA/kg'])][SalidA1],
                                             'Cultivation':diccionarios_colores[len(variables_items_ordened['Cultivation'])][SalidA2],
                                             'Coffee Variety':diccionarios_colores[len(variables_items_ordened['Coffee Variety'])][SalidA3],
                                             'Genomic DNA kit':diccionarios_colores[len(variables_items_ordened['Genomic DNA kit'])][SalidA4],
                                             'Drying Time (Days)':diccionarios_colores[len(variables_items_ordened['Drying Time (Days)'])][SalidA5],
                                             'Postharvest Processing':diccionarios_colores[len(variables_items_ordened['Postharvest Processing'])][SalidA6]}
            #print(asignaciones_colors_variables)  # claves con una lista de colores, cada variable tiene sus colores 
            
            var_asig_colors = {}
            for i in variables:
                var_children = metadata[i].unique().tolist()
                #---
                try:
                    var_children.sort(key=float)
                except ValueError:
                    var_children.sort(key=str)
                var_asig_colors[i] = dict(zip(var_children, asignaciones_colors_variables[i]))
                
            
            #print(var_asig_colors) # cada clave con su color asignado y ordenados, usar este para
            #-------------------------------------------------------
            limite = 100
            if tipo_kits_1 == 'Both kits':
                inters = 0.7
                AnchO = 9
                if width_linea == 0:
                    dendogramaX2 = 0.0000001
                if width_linea > 0:
                    dendogramaX2 = ancho_dendo # 0.17
                stackedX2 = 0.5 # ancho stackplot
                aumento = 0.002
                cuadro = 0.5
                dendogramaY2 = 1
                ancho_x2 = 0.5
                #----------------------------
                bbox_X, bbox_Y = despazamiento[displacement_leye], 1.053 # bbox_X, bbox_Y = 1.45, 1.055
                x_pos_leyendas, y_pos_leyendas = -1.04, 0
                y_pos_leyendas_intervalo = 0.035
                #------    sizes   ----------
                xlabel = sample_text+1
                label_25_75_100 = sample_text
                markersize1 = AnchO-5.2
                texttitle_size = AnchO-1
                label_legend1 = AnchO-2.5
                text_variables = varia_text # AnchO-2.6 
                text_samples = sample_text #AnchO-1.5
                label_legend1 = AnchO - 2.7
                markersize2 = AnchO - 4.1
                label_legend2 = AnchO - 2.5


            if tipo_kits_1 == 'DNPowerSoil':
                inters = 0.7
                inters = inters * (1.2/inters)
                AnchO = 10
                if width_linea == 0:
                    dendogramaX2 = 0.0000001
                if width_linea > 0:
                    dendogramaX2 = ancho_dendo # 0.17
                stackedX2 = (0.8 * len(DNMicrobial)) / len(name_sample) # ancho stackplot
                aumento = 0.002
                cuadro = 0.5
                dendogramaY2 = (1 * len(DNMicrobial)) /  len(name_sample)
                ancho_x2 = 0.25
                #----------------------------
                bbox_X, bbox_Y = despazamiento[displacement_leye] * 1.07, 1.1 # 1.565, 1.11
                x_pos_leyendas, y_pos_leyendas = -1.67,  -0.02
                y_pos_leyendas_intervalo = 0.065
                #------    sizes   ----------
                xlabel = sample_text+1
                label_25_75_100 = sample_text
                markersize1 = AnchO-5.2
                texttitle_size = AnchO-1
                label_legend1 = AnchO-2.5
                text_variables = varia_text + 1 # AnchO-2.6 
                text_samples = sample_text + 1 # AnchO-1.5
                label_legend1 = AnchO - 2.7
                markersize2 = AnchO - 4.1
                label_legend2 = AnchO - 2.5

            if tipo_kits_1 == 'DNMicrobial':
                inters = 0.7
                inters = inters * (1.2/inters)
                AnchO = 10
                if width_linea == 0:
                    dendogramaX2 = 0.0000001
                if width_linea > 0:
                    dendogramaX2 = ancho_dendo # 0.17
                stackedX2 = (0.8 * len(DNMicrobial)) / len(name_sample) # ancho stackplot
                aumento = 0.002
                cuadro = 0.5
                dendogramaY2 = (1 * len(DNMicrobial)) /  len(name_sample)
                ancho_x2 = 0.25
                #----------------------------
                bbox_X, bbox_Y = despazamiento[displacement_leye] * 1.07, 1.1# 1.565, 1.11
                x_pos_leyendas, y_pos_leyendas = -1.67,  -0.02
                y_pos_leyendas_intervalo = 0.065
                #------    sizes   ----------
                xlabel = sample_text+1
                label_25_75_100 = sample_text
                markersize1 = AnchO-5.2
                texttitle_size = AnchO-1
                label_legend1 = AnchO-2.5
                text_variables = varia_text + 1 # AnchO-2.6 
                text_samples = sample_text + 1 # AnchO-1.5
                label_legend1 = AnchO - 2.7
                markersize2 = AnchO - 4.1
                label_legend2 = AnchO - 2.5
            
            
            #-------------------------------------------------------
            inters_sum = (len(ivl) + 1) * inters
            bar_width = (limite - inters_sum) / len(ivl)
            radio = bar_width/2
            centros_circles = {}
            pos_inters = inters
            for i in ivl:
                centros_circles[i] = pos_inters + radio
                pos_inters += (inters + bar_width)
            # ajuste de los ejes para que concuerde con el stack plot
            inicioY = -abs(ivw * (list(centros_circles.values())[0] / limite) - np.min(np.array(ylines)))
            final_Y = abs(ivw * (list(centros_circles.values())[-1] / limite) - np.max(np.array(ylines))) + ivw



            #-----------------------------------------------------
            #-----------------------------------------------------
            #-----------------------------------------------------
            #-----------------------------------------------------

            mpl.rcParams.update(mpl.rcParamsDefault)
            fig = plt.figure(figsize = (AnchO, AnchO * cuadro))

            ax0 = fig.add_axes([0, 0, dendogramaX2, dendogramaY2])
            ax0.set_facecolor('none')
            ax0.set_xlim([dvw, 0])
            ax0.set_ylim([inicioY, final_Y])

            ax0.set_yticks(iv_ticks)

            ax0.yaxis.set_ticks_position('right')

            for line in ax0.get_yticklines():
                line.set_visible(False)

            colors_used = _remove_dups(color_list)
            color_to_lines = {}
            for color in colors_used:
                color_to_lines[color] = []
            for (xline, yline, color) in zip(xlines, ylines, color_list):
                color_to_lines[color].append(list(zip(xline, yline)))

            colors_to_collections = {}
            for color in colors_used:
                coll = matplotlib.collections.LineCollection(color_to_lines[color], colors=(color,))
                colors_to_collections[color] = coll

            for color in colors_used:
                if color != above_threshold_color:
                    ax0.add_collection(colors_to_collections[color])

            if above_threshold_color in colors_to_collections:
                    ax0.add_collection(colors_to_collections[above_threshold_color])

            for e, i in enumerate(ax0.collections):
                i.set_color(Set1[e])
                
            ax0.yaxis.set_visible(False)

            ax0.axis('off')
            
            #-----------------------------------------------------
            #-----------------------------------------------------
            #-----------------------------------------------------
            #-----------------------------------------------------

            ax1 = fig.add_axes([dendogramaX2 + aumento, 0, stackedX2, dendogramaY2])
            ax1.set_facecolor('none')
            ax1.set_xlim([0, 100])
            ax1.set_ylim([0, limite])

            for cc in mues:# x1     y1        x2      y2
                n = 0
                for i, t in zip(list(mues[cc].values()), list(mues[cc].keys())):
                    ax1.add_patch(Rectangle((n, centros_circles[cc]-radio), i, bar_width, fc =colors2[t],  ec ='white', lw = 0,
                                            url = 'https://www.google.com/search?client=firefox-b-d&q='+t.lower()+'+coffee'))
                    n += i


            ax1.set_xlabel('') # ('Relative abundance (%)', fontsize= xlabel, fontname="Open Sans", weight = 'bold')
            ax1.set_yticklabels([])

            for d in [25, 50, 75]:
                Lista = list(range(0, 100))
                ax1.plot(np.repeat(d, 100+1),np.array(Lista+[Lista[-1]+3]), 
                         marker='o', markeredgewidth=0, zorder=0, linestyle='-',
                                     markersize=0, color='black', linewidth=0.3, alpha = 1)
            plt.xticks([25, 50, 75, 100], ['25', '50', '75', '100'], size= label_25_75_100, rotation=0, ha = 'right')

            ax1.xaxis.set_ticks_position('top')

            plt.gca().tick_params(bottom=False, right=False, top=True, left=False, width = 0.3, length=3, color='black')

            plt.gca().spines['left'].set_color(None)
            plt.gca().spines['bottom'].set_color(None)
            plt.gca().spines['right'].set_color(None)
            plt.gca().spines['top'].set_color(None)

            ###  legend

            proxies = []
            for cat in colors2:
                proxy = mpl.lines.Line2D([0], [0], linestyle='',
                                         c=colors2[cat], marker='s', alpha = 1,
                                         markersize = markersize1,
                                         markeredgecolor=colors2[cat], linewidth = 0)
                proxies.append(proxy)

            ax1.legend(proxies, list(colors2.keys()), title= "$\\bf{"+tax+"}$", title_fontsize = texttitle_size, numpoints=1, loc=2, ncol = num_cols,
                      handletextpad=0.7, handlelength = 0.3, labelspacing = 0.5, columnspacing = 1,
                               borderpad = 0.5, edgecolor="none", prop={'style':'italic', 'size': label_legend1},
                              bbox_to_anchor=(bbox_X, bbox_Y))

            ax1.xaxis.set_ticks_position("top")
            ax1.set_title('Relative abundance (%)', fontsize= xlabel, fontname="Open Sans", weight = 'bold')
            
            
            ###

            #-----------------------------------------------------
            #-----------------------------------------------------
            #-----------------------------------------------------
            #-----------------------------------------------------
            
            hide_variables = list(set(variables) - set(list(Variables_metadatA))) # lista con variables que no se quieren mostrar
            
            
            
            if 'None' in Variables_metadatA:
                with open('plots_asv/taxonomy/dict_variable_element_colors.json', 'w') as fp:
                    json.dump({'None':'None'}, fp)
            else:
                DDD = {}
                for v in Variables_metadatA:
                    DDD[v] = dict(zip(variables_items_ordened[v], asignaciones_colors_variables[v]))
                    
                #print(DDD) # diccionario ajustable a la  cantidad de variables elegidas
            
                with open('plots_asv/taxonomy/dict_variable_element_colors.json', 'w') as fp:
                    json.dump(DDD, fp)
                
            VVVVV = {}        
            for v in variables:
                VVVVV[v] = dict(zip(variables_items_ordened[v], asignaciones_colors_variables[v]))
            with open('plots_asv/taxonomy/dict_variable_element_colors_ALL.json', 'w') as fp:
                json.dump(VVVVV, fp)
             
                    
                    
                    
                    
            ax2 = fig.add_axes([dendogramaX2 + stackedX2 + aumento, 0, ancho_x2, dendogramaY2])
            ax2.set_facecolor('none')
            ax2.set_xlim([0, 100])
            ax2.set_ylim([0, 100])

            pos_variable_title = {}
            num = radio*2
            for i in variables:
                if i in hide_variables:
                    pass
                else:
                    for e, cc in enumerate(mues):
                        ax2.add_patch(Circle((num, centros_circles[cc]), radio, fc = var_asig_colors[i][correspondencia_sam_vars[cc][i]],  ec ='white', lw = 0, zorder = 2))
                    pos_variable_title[i] = num
                    num += (radio*2.2)

            for VaR in pos_variable_title:
                ax2.text(pos_variable_title[VaR], 100, ' '+VaR, fontsize = text_variables, color = 'black', ha='center', va = 'bottom', rotation=90, fontfamily = family_axis) 


            for sam in centros_circles:
                #ax2.text(num - radio, centros_circles[sam], ' '+sam, fontsize = text_samples, color = 'black', ha='left', va = 'center', rotation=0, fontfamily = family_axis)
                
                
                if cambiar_sample == 'Code1':
                    ax2.annotate(' '+sam, xy = (num - radio, centros_circles[sam]), xytext = (num - radio, centros_circles[sam]), xycoords = 'data', textcoords = 'data',
                                 url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[sam],
                                    color = 'black', size = text_samples, ha='left', va = "center", fontfamily = family_axis,
                            bbox = dict(boxstyle="Square", pad = 0.1, alpha = 0.1, fc = 'whitesmoke', ec = "none", url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[sam]))
                if cambiar_sample == 'Code2':
                    ax2.annotate(' '+name_code[sam], xy = (num - radio, centros_circles[sam]), xytext = (num - radio, centros_circles[sam]), xycoords = 'data', textcoords = 'data',
                                 url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[sam],
                                    color = 'black', size = text_samples, ha='left', va = "center", fontfamily = family_axis,
                            bbox = dict(boxstyle="Square", pad = 0.1, alpha = 0.1, fc = 'whitesmoke', ec = "none", url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[sam]))
                    
                if cambiar_sample == 'Code3':
                    ax2.annotate(' '+name_code2[sam], xy = (num - radio, centros_circles[sam]), xytext = (num - radio, centros_circles[sam]), xycoords = 'data', textcoords = 'data',
                                 url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[sam],
                                    color = 'black', size = text_samples, ha='left', va = "center", fontfamily = family_axis,
                            bbox = dict(boxstyle="Square", pad = 0.1, alpha = 0.1, fc = 'whitesmoke', ec = "none", url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[sam]))
                
                
            #

            ly = y_pos_leyendas
            leyendas = []
            for cat in pos_variable_title:
                labels = []
                proxies = []

                proxyy = mpl.lines.Line2D([-10], [0], linestyle='none',
                                             c='black', marker='', alpha = 1,
                                             markersize = 5, 
                                             markeredgecolor='black', linewidth = 0)
                proxies.append(proxyy)
                labels.append(' '.join(["$\\bf{"+i+"}$" for i in cat.split(' ')])+':') # separar las palabras para editar cada una

                for vari in var_asig_colors[cat]:
                    proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                             c=var_asig_colors[cat][vari], marker='o', alpha = 1,
                                             markersize = markersize2,
                                             markeredgecolor=var_asig_colors[cat][vari], linewidth = 0)
                    labels.append(vari)
                    proxies.append(proxy)


                leg = ax2.legend(proxies, labels, title_fontsize = 8, numpoints=1, loc=2, ncol = len(labels), facecolor = 'none',
                          handletextpad=0.7, handlelength = 0.3, labelspacing = 0, columnspacing = 1.5,
                                   borderpad = 0.13, edgecolor="none", prop={'size':label_legend2},
                                  bbox_to_anchor=(x_pos_leyendas, ly)) 

                leyendas.append(leg)
                ly -= y_pos_leyendas_intervalo



            for l in leyendas:
                ax2.add_artist(l)

            
            
            ax2.axis('off')

            plt.close()
            
            
            
            # --- save
            ancho_plot_save = str(71)
            png1 = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
            png1.style.button_color = 'gold'
            output1 = widgets.Output()
            def button_clicked1(b):
                with output1:
                    clear_output(True)
                    if filename_plot.value is '':
                        nombre_grafico = 'Taxonomy_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                    if filename_plot.value is not '':
                        nombre_grafico = filename_plot.value
                    fig.savefig('plots_asv/taxonomy/'+nombre_grafico+'.png', dpi = 900, bbox_inches= 'tight')
                    
            png1.on_click(button_clicked1)
            #----
            jpeg1 = widgets.Button(description="JPEG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
            jpeg1.style.button_color = 'gold'
            output3 = widgets.Output()
            def button_clicked3(b):
                with output3:
                    clear_output(True)
                    if filename_plot.value is '':
                        nombre_grafico = 'Taxonomy_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                    if filename_plot.value is not '':
                        nombre_grafico = filename_plot.value
                    fig.savefig('plots_asv/taxonomy/'+nombre_grafico+'.jpeg', dpi = 900, bbox_inches= 'tight')
                    
            jpeg1.on_click(button_clicked3)
            #----
            svg1 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
            svg1.style.button_color = 'gold'
            output2 = widgets.Output()
            def button_clicked2(b):
                with output2:
                    clear_output(True)
                    if filename_plot.value is '':
                        nombre_grafico = 'Taxonomy_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                    if filename_plot.value is not '':
                        nombre_grafico = filename_plot.value
                    fig.savefig('plots_asv/taxonomy/'+nombre_grafico+'.svg', dpi = 900, bbox_inches= 'tight')
                    
            svg1.on_click(button_clicked2)
            #----
            pdf1 = widgets.Button(description="PDF", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
            pdf1.style.button_color = 'gold'
            output4 = widgets.Output()
            def button_clicked4(b):
                with output4:
                    clear_output(True)
                    if filename_plot.value is '':
                        nombre_grafico = 'Taxonomy_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                    if filename_plot.value is not '':
                        nombre_grafico = filename_plot.value
                    fig.savefig('plots_asv/taxonomy/'+nombre_grafico+'.pdf', dpi = 900, bbox_inches= 'tight')
                    
            pdf1.on_click(button_clicked4)
            #----
            
            
            params = widgets.Button(description="SAVE", icon = 'fa-save', layout=Layout(width='70px'))
            params.style.button_color = 'cyan'
            output5 = widgets.Output()
            def button_clicked5(b):
                with output5:
                    clear_output(True)
                    
                    if filename.value is '':
                        nombre_archivo = 'Taxonomy_chart_params_ASVs.txt'
                    if filename.value is not '':
                        nombre_archivo = filename.value
                    
                    with open('plots_asv/taxonomy/'+nombre_archivo, 'w') as fq:
                        fq.write('#Saved parameters\n')
                        fq.write('Clustering:ASVs\n' )
                        """
                        diccionario para guardar en un archivo los parametros usados para volver a reproducirlos posteriormente
                        """
                        parametros = {'Kit:':tipo_kits_1,
                                      'Level:':tax,
                                      'Limit:':str(PercentagE),
                                      'Variables:':', '.join(list(Variables_metadatA)),
                                      'Method:':METODO,
                                      'Metric:':METRICA,
                                      'Show:':para_mostrar,
                                      'Chart size:':tam_plot1,
                                      'Line width:':width_linea,
                                      'Sample text:':sample_text,
                                      'Variable text:':varia_text,
                                      'Axis font:':family_axis,
                                      'Legend pos':displacement_leye,
                                      'Columns:':num_cols,
                                      'Dendrogram width:':ancho_dendo,
                                      'BAR COLORS:':', '.join(list(Multiple_colorS)),
                                      '\n###':' Variable colors in metadata ###\n',
                                      'Location:':SalidA0+' = ('+', '.join(diccionarios_colores[len(variables_items_ordened['Location'])][SalidA0])+')',
                                      'ug OTA/kg:':SalidA1+' = ('+', '.join(diccionarios_colores[len(variables_items_ordened['ug OTA/kg'])][SalidA1])+')',
                                      'Cultivation:':SalidA2+' = ('+', '.join(diccionarios_colores[len(variables_items_ordened['Cultivation'])][SalidA2])+')',
                                      'Coffee Variety:':SalidA3+' = ('+', '.join(diccionarios_colores[len(variables_items_ordened['Coffee Variety'])][SalidA3])+')',
                                      'Genomic DNA kit:':SalidA4+' = ('+', '.join(diccionarios_colores[len(variables_items_ordened['Genomic DNA kit'])][SalidA4])+')',
                                      'Drying Time (Days):':SalidA5+' = ('+', '.join(diccionarios_colores[len(variables_items_ordened['Drying Time (Days)'])][SalidA5])+')',
                                      'Postharvest Processing:':SalidA6+' = ('+', '.join(diccionarios_colores[len(variables_items_ordened['Postharvest Processing'])][SalidA6])+')'}
                        for w in parametros:
                            fq.write(w+str(parametros[w])+'\n')

            params.on_click(button_clicked5)
            #------------------------------------
            


            def set_locus():
                
                for e, i in enumerate(ax0.collections):
                    i.set_linewidth(width_linea)
                
                
                if tipo_kits_1 == 'Both kits':
                    fig.set_size_inches((tam_plot1-1, (tam_plot1-1) * 0.5))
                if tipo_kits_1 == 'DNPowerSoil':
                    fig.set_size_inches((tam_plot1, tam_plot1 * 0.5))
                if tipo_kits_1 == 'DNMicrobial':
                    fig.set_size_inches((tam_plot1, tam_plot1 * 0.5))
                
                
                
                display(fig)

            
            #------------------------------------
            
            Sample_Select = widgets.Dropdown(options=list(KITS[tipo_kits_1]), disabled=False,
                                     layout=Layout(width='160px', height='27px'))

            blanca33 = widgets.Button(layout=Layout(width='10px', height='2px'), disabled=True)
            blanca33.style.button_color = 'white'

            explore_title = widgets.HTML('<font color = grey> <i class="fa fa-eye fa-2x fa-fw"></i> <b style="font-size:0.8vw">EXPLORE SAMPLE INFORMATION</b>')

            def samsel(Sample_Select):
                bot0 = widgets.Button(description=Sample_Select, icon = 'fa-hand-pointer-o', layout=Layout(width='160px', height='27px'))
                bot0.style.button_color = 'lightgreen'
                bot0.style.font_weight = 'bold'
                outbot0 = widgets.Output()

                display(bot0, outbot0)
                def button_bot0(b):
                    with outbot0:
                        clear_output(True)
                        SELECTSAM(SAM_SELECT = Sample_Select)

                bot0.on_click(button_bot0)
            samselOUT = widgets.interactive_output(samsel, {'Sample_Select':Sample_Select})

            exploracion_de_muestras = VBox([explore_title, HBox([VBox([blanca33, Sample_Select]), samselOUT, HBox([VBox([blanca33, HBox([blanca2, widgets.Label('Intersections:'), Chord])])])])])

            exploracion_de_muestras_box = Box(children=[exploracion_de_muestras], layout=Layout(border='1px solid gainsboro', width='580px', height='80px'))
            
            
            
            items_save_1 = VBox([HBox([widgets.HTML('<font color = grey> <b style="font-size:0.7vw">SAVE CHART: </b>'), blanca, filename_plot]),
                                 HBox([blanca, widgets.Label('Formats:'), png1, jpeg1, svg1, pdf1]),
                         HBox([blanca, widgets.Label('Chart parameters:'), filename, params])])
            items_save_1_box = Box(children=[items_save_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_save_1.children) * 34))+'px'))
            
            display(HBox([VBox([items_save_1_box, cuenta_colors_tax]), exploracion_de_muestras_box]))
            
            display(Chord_out)
            
            set_locus()
            
            #uno.value = '<i class="fa fa-spinner fa-2x fa-fw"></i> <font color = black>  <h style="font-size:0.6vw"> Finishied.</font>'
            uno.value = '<i class="fa fa-spinner fa-2x fa-fw"></i>'

        OUT_TaxonomY = widgets.interactive_output(TaxonomY, {'tipo_kits_1':tipo_kits_1, 'tax':tax, 'PercentagE':PercentagE,
                                                             'tam_plot1':tam_plot1, 'width_linea':width_linea, 'METODO':METODO, 'METRICA':METRICA,
                                                             'Variables_metadatA':Variables_metadatA, 'Multiple_colorS':Multiple_colorS, 'para_mostrar':para_mostrar,
                                                             'family_axis':family_axis, 'displacement_leye':displacement_leye, 'num_cols':num_cols, 'SalidA0':SalidA0,
                                                             'SalidA1':SalidA1, 'SalidA2':SalidA2, 'SalidA3':SalidA3, 'SalidA4':SalidA4, 'SalidA5':SalidA5, 'SalidA6':SalidA6,
                                                             'varia_text':varia_text, 'sample_text':sample_text, 'ancho_dendo':ancho_dendo, 'cambiar_sample':cambiar_sample}) 
        
        items_plot = VBox([HBox([blanca, widgets.Label('Method:'), METODO]),
                                HBox([blanca, widgets.Label('Metric:'), METRICA]),
                                HBox([blanca, widgets.Label('Show:'), para_mostrar]),
                              tam_plot1,
                                width_linea,
                           sample_text,
                           varia_text,
                               HBox([blanca, widgets.Label('Axis font:'), family_axis]),
                               displacement_leye,
                               HBox([blanca, widgets.Label('Columns:'), num_cols]),
                           HBox([widgets.Label('Rename sample:'), cambiar_sample]),
                          ancho_dendo])
        items_plot_box = Box(children=[items_plot], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_plot.children) * 31))+'px'))
        
        
        taxtaxtax = HBox([VBox([widgets.HTML('<font color = #1976d2> <b style="font-size:0.8vw">KITS</b>'), tipo_kit_1_box,
                              widgets.HTML('<font color = grey> <b style="font-size:0.8vw">LEVEL</b>'), HBox([VBox([tipo_tax_box]), blanca2, PercentagE_box]),
                                widgets.HTML('<font color = grey> <b style="font-size:0.8vw">VARIABLES</b>'), Variables_metadatA_box]),
                          VBox([widgets.HTML('<font color = grey> <i class="fa fa-cog fa-2x fa-fw"></i> <b style="font-size:0.8vw">PLOT SETTINGS</b>'),
                               
                               items_plot_box,
                                HBox([VBox([widgets.HTML('<font color = grey> <b style="font-size:0.6vw">BAR COLORS</b>'), Multiple_colorS, ayuda3]),
                               VBox([widgets.HTML('<font color = grey> <b style="font-size:0.6vw">VARIABLE COLORS</b>'), Colores_para_variableS])])]),
                                
                                
                          VBox([OUT_TaxonomY])])
        
        
        display(taxtaxtax)


TAXONOMY_button.on_click(button_clicked)





TAXONOMY_ANALYSIS = VBox([HBox([TAXONOMY_button, uno]), TAXONOMY_output])














# # BETA DIVERSIRY: NO PHYLOGENETIC
# ## BRAY CURTIS DISTANCES ANALYSIS AND THEN HIERARCHICAL CLUSTERING USING UPGMA (from scipy.cluster.hierarchy import average)  




from numpy import zeros
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import squareform

def bray_curtis_distance(table, sample1_id, sample2_id):
    numerator = 0
    denominator = 0
    sample1_counts = table[sample1_id]
    sample2_counts = table[sample2_id]
    for sample1_count, sample2_count in zip(sample1_counts, sample2_counts):
        numerator += abs(sample1_count - sample2_count)
        denominator += sample1_count + sample2_count
    return numerator / denominator
def table_to_distances(table, pairwise_distance_fn):
    sample_ids = table.columns
    num_samples = len(sample_ids)
    data = zeros((num_samples, num_samples))
    for i, sample1_id in enumerate(sample_ids):
        for j, sample2_id in enumerate(sample_ids[:i]):
            data[i,j] = data[j,i] = pairwise_distance_fn(table, sample1_id, sample2_id)
    return data, list(sample_ids)

###############  obtenido a partir de SKBIO
from functools import partial
def _preprocess_input(distance_matrix, grouping, column):
    if isinstance(grouping, pd.DataFrame):
        if column is None:
            raise ValueError(
                "Must provide a column name if supplying a DataFrame.")
        else:
            grouping = _df_to_vector(distance_matrix, grouping, column)
    elif column is not None:
        raise ValueError(
            "Must provide a DataFrame if supplying a column name.")

    sample_size = distance_matrix.shape[0]
    if len(grouping) != sample_size:
        raise ValueError(
            "Grouping vector size must match the number of IDs in the "
            "distance matrix.")

    # Find the group labels and convert grouping to an integer vector
    # (factor).
    groups, grouping = np.unique(grouping, return_inverse=True)
    num_groups = len(groups)

    if num_groups == len(grouping):
        raise ValueError(
            "All values in the grouping vector are unique. This method cannot "
            "operate on a grouping vector with only unique values (e.g., "
            "there are no 'within' distances because each group of objects "
            "contains only a single object).")
    if num_groups == 1:
        raise ValueError(
            "All values in the grouping vector are the same. This method "
            "cannot operate on a grouping vector with only a single group of "
            "objects (e.g., there are no 'between' distances because there is "
            "only a single group).")

    tri_idxs = np.triu_indices(sample_size, k=1)
    distances = squareform(distance_matrix, force='tovector', checks=False) # CAMBIOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    return sample_size, num_groups, grouping, tri_idxs, distances

def _df_to_vector(distance_matrix, df, column):

    if column not in df:
        raise ValueError("Column '%s' not in DataFrame." % column)

    grouping = df.reindex(df.index, axis=0).loc[:, column] # CAMBIOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    if grouping.isnull().any():
        raise ValueError(
            "One or more IDs in the distance matrix are not in the data "
            "frame.")
    return grouping.tolist()

def _compute_f_stat(sample_size, num_groups, tri_idxs, distances, group_sizes,
                    s_T, grouping):
    """Compute PERMANOVA pseudo-F statistic."""
    # Create a matrix where objects in the same group are marked with the group
    # index (e.g. 0, 1, 2, etc.). objects that are not in the same group are
    # marked with -1.
    grouping_matrix = -1 * np.ones((sample_size, sample_size), dtype=int)
    for group_idx in range(num_groups):
        within_indices = _index_combinations(
            np.where(grouping == group_idx)[0])
        grouping_matrix[within_indices] = group_idx

    # Extract upper triangle (in same order as distances were extracted
    # from full distance matrix).
    grouping_tri = grouping_matrix[tri_idxs]

    # Calculate s_W for each group, accounting for different group sizes.
    s_W = 0
    for i in range(num_groups):
        s_W += (distances[grouping_tri == i] ** 2).sum() / group_sizes[i]

    s_A = s_T - s_W
    return (s_A / (num_groups - 1)) / (s_W / (sample_size - num_groups))

def _index_combinations(indices):
    # Modified from http://stackoverflow.com/a/11144716
    return np.tile(indices, len(indices)), np.repeat(indices, len(indices))
def _run_monte_carlo_stats(test_stat_function, grouping, permutations):
    """Run stat test and compute significance with Monte Carlo permutations."""
    if permutations < 0:
        raise ValueError(
            "Number of permutations must be greater than or equal to zero.")

    stat = test_stat_function(grouping)

    p_value = np.nan
    if permutations > 0:
        perm_stats = np.empty(permutations, dtype=np.float64)

        for i in range(permutations):
            perm_grouping = np.random.permutation(grouping)
            perm_stats[i] = test_stat_function(perm_grouping)

        p_value = ((perm_stats >= stat).sum() + 1) / (permutations + 1)

    return stat, p_value
def _build_results(method_name, test_stat_name, sample_size, num_groups, stat,
                   p_value, permutations):
    """Return ``pandas.Series`` containing results of statistical test."""
    return pd.Series(
        data=[method_name, test_stat_name, sample_size, num_groups, stat,
              p_value, permutations],
        index=['method name', 'test statistic name', 'sample size',
               'number of groups', 'test statistic', 'p-value',
               'number of permutations'],
        name='%s results' % method_name)


def permanova(distance_matrix, grouping, column=None, permutations=999):
    sample_size, num_groups, grouping, tri_idxs, distances = _preprocess_input(
        distance_matrix, grouping, column)

    # Calculate number of objects in each group.
    group_sizes = np.bincount(grouping)
    s_T = (distances ** 2).sum() / sample_size

    test_stat_function = partial(_compute_f_stat, sample_size, num_groups,
                                 tri_idxs, distances, group_sizes, s_T)
    stat, p_value = _run_monte_carlo_stats(test_stat_function, grouping,
                                           permutations)

    return _build_results('PERMANOVA', 'pseudo-F', sample_size, num_groups,
                          stat, p_value, permutations)


from matplotlib.patches import Ellipse, Circle, Rectangle





import plotly.graph_objects as go
import plotly.io





# dos niveles
import datetime

progreso1_no_fil = widgets.HTML('<font color = black> <i class="fa fa-spinner fa-2x fa-fw"></i> </font>')
progreso2_no_fil = widgets.HTML()
estatico_no_fil = widgets.HTML('<font color = limegreen> <b style="font-size:0.5vw">Processed samples : </b>')
estatico_no_fil = HBox([progreso1_no_fil, estatico_no_fil, progreso2_no_fil])
estatico_no_fil_box = Box(children=[estatico_no_fil], layout=Layout(border='1px solid limegreen', width='585px', height='32px'))

BETA_button = widgets.Button(description="PROCESS AND VISUALIZE", icon = 'fa-eye', layout=Layout(width='590px'))
BETA_button.style.button_color = 'gainsboro' #'deepskyblue'
BETA_button.style.font_weight = 'bold'
BETA_output = widgets.Output()

def button_clicked(b):
    import time
    from scipy.cluster.hierarchy import dendrogram, average
    
    
    NCBI_RDP_SILVA_SUMMARY = pd.read_csv('tablas/ASVs_NCBI_RDP_SILVA_SUMMARY.txt', sep = '\t')
    
    rarefaction_collapse_specie = pd.pivot_table(NCBI_RDP_SILVA_SUMMARY[['Species'] + name_sample], values = name_sample, index = ['Species'], aggfunc = sum).reset_index()
    
    with open('plots_asv/taxonomy/dict_variable_element_colors_ALL.json', 'r') as fp:
        dict_variable_element_colors = json.load(fp)

    asv_species_dict = {}
    asv_species_lista = []
    for i in NCBI_RDP_SILVA_SUMMARY.Species.drop_duplicates():
        w = NCBI_RDP_SILVA_SUMMARY[NCBI_RDP_SILVA_SUMMARY.Species == i]
        sumas = np.sum(w.iloc[:, 3:].values, axis = 1)
        #sumas = np.sum(w.values[0].tolist()[3:])
        w['sum'] = sumas
        w = w.sort_values(by =['sum'],ascending=False).reset_index(drop=True)
        ef = w['Entry'].tolist()[0]
        asv_species_dict[ef] = i
        asv_species_lista.append([ef, i])
    unicos = DataFrame(asv_species_lista, columns = ['Entry', 'Species'])
    
    MERGE_UNICOS_SPECIES = unicos.merge(rarefaction_collapse_specie, on = 'Species', how = 'left')
    
    with BETA_output:
        clear_output(True)
        
        progreso1_no_fil.value = '<font color = black> <i class="fa fa-spinner fa-pulse fa-2x fa-fw"></i> </font>'
        
        
        PCA_KITS = {}
        for k in KITS:
            if k == 'Both kits':
                NEW_rarefaction_collapse_specie = MERGE_UNICOS_SPECIES[['Species'] + KITS['Both kits']]
            if k == 'DNPowerSoil':
                NEW_rarefaction_collapse_specie = MERGE_UNICOS_SPECIES[['Species'] + KITS['DNPowerSoil']]
            if k == 'DNMicrobial':
                NEW_rarefaction_collapse_specie = MERGE_UNICOS_SPECIES[['Species'] + KITS['DNMicrobial']]

            ONE = NEW_rarefaction_collapse_specie.set_index('Species')
            
            bc_matrix, id_samples = table_to_distances(ONE, bray_curtis_distance)
            
            
            SAMple = DataFrame(id_samples, columns = ['Name Sample']).merge(metadata, on = 'Name Sample', how = 'left')
            SAMple =SAMple.set_index('Name Sample')
            DfDf = DataFrame(bc_matrix, columns = id_samples)
            DfDf.insert(loc = 0, column='Name Sample', value= id_samples)
            DfDf = DfDf.merge(metadata, on = 'Name Sample', how = 'left')
            DfDf = DfDf.set_index('Name Sample')

            SS = StandardScaler()
            DfDf[id_samples] = SS.fit_transform(DfDf[id_samples].astype('float'))

            pca3 = PCA(n_components=3)
            pca_3 = pca3.fit_transform(DfDf[id_samples].astype('float'))

            pca_kits_variables = {}

            for variable in variables:
                progreso2_no_fil.value = '<font color = red> <b style="font-size:0.5vw">'+k+' = '+variable+'</b>'
                
                if len(DfDf[variable].unique()) == 1:
                    pass
                else:
                    DF_3 = pd.DataFrame({'PCA1':pca_3[:,0], 'PCA2':pca_3[:, 1], 'PCA3':pca_3[:, 2], 'clase':DfDf[variable]})
                    Permanova = permanova(bc_matrix, SAMple, variable)['p-value']

                    bc_matrix_condensed_ASVs = squareform(bc_matrix, force='tovector', checks=False)
                    ZZZ = average(bc_matrix_condensed_ASVs)

                    dendro = dendrogram(ZZZ,
                                orientation='left',  no_plot = True,
                                labels=id_samples,
                                distance_sort='descending',
                                        leaf_font_size = 100,
                                show_leaf_counts=True)
                    #-----------------------------------------
                    ivl = dendro['ivl']
                    color_list = dendro['color_list']
                    Z = np.asarray(dendro['dcoord'], order='c')
                    mh = max(Z[:, 2])
                    above_threshold_color = 'b'
                    # Independent variable plot width
                    ivw = len(ivl) * 10
                    # Dependent variable plot height
                    dvw = mh + mh * 0.05
                    iv_ticks = np.arange(5, len(ivl) * 10 + 5, 10)
                    xlines = dendro['dcoord']
                    ylines = dendro['icoord']

                    category = dict_variable_element_colors[variable]
                    col = DF_3['clase'].map(category)


                pca_kits_variables[variable] = {'df_3':DF_3, 'permanova':Permanova, 'ivl':ivl, 'ivw':ivw, 'dvw':dvw, 'category':category,
                                                'iv_ticks':iv_ticks, 'xlines':xlines, 'ylines':ylines, 'color_list':color_list, 'above_threshold_color':above_threshold_color}
                time.sleep(0.1)
                
                
            PCA_KITS[k] = pca_kits_variables
            
        del DF_3, Permanova, ivl, ivw, dvw, category, iv_ticks, xlines, ylines, color_list, above_threshold_color, id_samples
        
        clear_output(True)
        progreso1_no_fil.value = '<font color = black> <i class="fa fa-spinner fa-2x fa-fw"></i> </font>'
        progreso2_no_fil.value = '<font color = red> <b style="font-size:0.5vw">'+k+' = '+variable+'</b>'
        
        
        mpl.rcParams.update(mpl.rcParamsDefault)

        fig = plt.figure(figsize=(6, 4))

        ax = fig.add_axes([0, 0, 1, 1])

        plt.gca().tick_params(which='major', width = 2, length=4, color='gainsboro')
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_color('gainsboro')
        plt.gca().spines['bottom'].set_color('gainsboro')
        plt.gca().spines['right'].set_color(None)
        plt.gca().spines['top'].set_color(None)

        plt.gca().set_xlabel('Contigs', fontsize=12, fontname='Open Sans', weight = 'bold')
        plt.gca().set_ylabel('Species', fontsize=12, fontname='Open Sans', weight = 'bold')
        plt.close()
        
        
        blanca = widgets.Button(layout=Layout(width='10px', height='25px'), disabled=True)
        blanca.style.button_color = 'white'
        blanca2 = widgets.Button(layout=Layout(width='1px', height='25px'), disabled=True)
        blanca2.style.button_color = 'white'
        
        
        tipo_kits_1 = widgets.ToggleButtons(options= ['Both kits'] + metadata[VARIABLE_KIT].unique().tolist(), value = 'Both kits', button_style = 'primary')
        tipo_kits_1.style.button_width = '170px'
        tipo_kits_1.style.font_weight = 'bold'
        tipo_kit_1_box = Box(children=[VBox([tipo_kits_1])], layout=Layout(border='1px solid #1976d2', width='180px', height='95px'))

        VariablE = widgets.ToggleButtons(options= variables, value = variables[0], button_style = 'warning')
        VariablE.style.button_width = '170px'
        tipo_variable_box = Box(children=[VBox([VariablE])], layout=Layout(border='1px solid #ff9800', width='180px', height='210px'))

        select_PCA = widgets.ToggleButtons(options=['2D', '3D', '3Di'], value = '2D')
        select_PCA.style.button_width = '50px'
        boton_select_PCA = Box(children=[select_PCA], layout= Layout(border='1px solid pink', width='170px', height='35px'))
        
        markers_point = widgets.ToggleButtons(options=['⏺', '⏹'], value = '⏺')
        markers_point.style.button_width = '31px'
        
        size_point = widgets.SelectionSlider(options=np.round(np.arange(0, 1.505, 0.05), 2), value=0.5,disabled=False,
                                              description = 'Marker size:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        
        alfa_point = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.7,disabled=False,
                                              description = 'Marker alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        alfa_linea = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.3,disabled=False,
                                              description = 'Line alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px')) 
        
        ver_num_muestra = widgets.ToggleButtons(options=['True', 'False'], value = 'False')
        ver_num_muestra.style.button_width = '55px'
        
        
        ver_dendograma = widgets.ToggleButtons(options=['True', 'False'], value = 'True')
        ver_dendograma.style.button_width = '55px'
        
        ver_var_name = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.1,disabled=False,
                                              description = 'Text alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        
        family = sorted(['Liberation Serif','Microsoft Sans Serif','Open Sans','Times New Roman','3ds Light','Calibri','Comic Sans MS',
                  'Arial','Courier New','Microsoft Yi Baiti','Lucida Console'])
        family_den = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                   layout = Layout(width='290px', height='25px'))

        ### ejes
        label_X = widgets.SelectionSlider(options=range(5, 31),value=12,disabled=False,
                                              description = 'X label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_X = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'X tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))

        label_Y = widgets.SelectionSlider(options=range(5, 31),value=12,disabled=False,
                                              description = 'Y label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_Y = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'Y tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))

        family_axis = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                   layout = Layout(width='290px', height='25px'))
        
        label_Z = widgets.SelectionSlider(options=range(5, 31),value=12,disabled=False,
                                              description = 'Z label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        
        ticklabels_Z = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'Z tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        
        
        filename = widgets.Text(value='', placeholder='File name.txt', description='', disabled=False, layout = Layout(width='195px', height='25px'))
        
        filename_plot = widgets.Text(value='', placeholder='Chart name', description='', disabled=False, layout = Layout(width='270px', height='25px'))
        
        cambiar_sample = widgets.ToggleButtons(options=['Code1', 'Code2', 'Code3'], value = 'Code1')
        cambiar_sample.style.button_width = '70px'

        
        def BetA(tipo_kits_1, VariablE, select_PCA):
            ElementoS = list(PCA_KITS[tipo_kits_1][VariablE]['category'].keys())
            
            DF_3 = PCA_KITS[tipo_kits_1][VariablE]['df_3']
            sam_clase = dict(zip(DF_3.index, DF_3.clase))
            Permanova = PCA_KITS[tipo_kits_1][VariablE]['permanova']
            
            ######

            VAR0 = 0
            VAR1 = 1
            VAR2 = 2
            VAR3 = 3
            VAR4 = 4
            VAR5 = 5
            VAR6 = 6

            width_cuadro = str(100)
            
            blanca3 = widgets.Button(layout=Layout(width='10px', height='2px'), disabled=True)
            blanca3.style.button_color = 'white'
            
            if VariablE == 'Location':
                for unoo in [variables[VAR0]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
            if VariablE == 'ug OTA/kg':
                for unoo in [variables[VAR1]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Cultivation':
                for unoo in [variables[VAR2]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Coffee Variety':
                for unoo in [variables[VAR3]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Genomic DNA kit':
                for unoo in [variables[VAR4]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Drying Time (Days)':
                for unoo in [variables[VAR5]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Postharvest Processing':
                for unoo in [variables[VAR6]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
            
            
            
            ivl = PCA_KITS[tipo_kits_1][VariablE]['ivl']
            ivw = PCA_KITS[tipo_kits_1][VariablE]['ivw']
            dvw = PCA_KITS[tipo_kits_1][VariablE]['dvw']
            iv_ticks = PCA_KITS[tipo_kits_1][VariablE]['iv_ticks']
            xlines = PCA_KITS[tipo_kits_1][VariablE]['xlines']
            ylines = PCA_KITS[tipo_kits_1][VariablE]['ylines']
            color_list = PCA_KITS[tipo_kits_1][VariablE]['color_list']
            above_threshold_color = PCA_KITS[tipo_kits_1][VariablE]['above_threshold_color']
            id_samples = KITS['Both kits']
            

            
            if tipo_kits_1 in ['DNPowerSoil', 'DNMicrobial'] and VariablE == 'Genomic DNA kit':
                print('!!! All values in the grouping vector are the same. !!!')
                
            else:
            
                if select_PCA in ('2D', '3D'):
                    import datetime
                    
                    def BetA2(markers_point, size_point, alfa_point, alfa_linea, ver_num_muestra, ver_dendograma, ver_var_name,
                              family_den, label_X, ticklabels_X, label_Y, ticklabels_Y, family_axis, label_Z, ticklabels_Z, cambiar_sample, SalidA0):
                        
                        category = dict(zip(ElementoS, diccionarios_colores[len(variables_items_ordened[VariablE])][SalidA0]))
                        col = DF_3['clase'].map(category)    
                        
                        if select_PCA == '2D':

                            LETRA_2pca = 10
                            plot_alto_pca2 = 6

                            mpl.rcParams.update(mpl.rcParamsDefault)

                            fig = plt.figure(figsize = (plot_alto_pca2, plot_alto_pca2))

                            ax = fig.add_axes([0, 0, 1, 0.8])
                            ax.set_aspect('equal')

                            ax.set_xlim(DF_3.PCA1.min() + (DF_3.PCA1.min() * 0.2), DF_3.PCA1.max() + (DF_3.PCA1.max() * 0.2))
                            ax.set_ylim(DF_3.PCA2.min() + (DF_3.PCA2.min() * 0.2), DF_3.PCA2.max() + (DF_3.PCA2.max() * 0.2))


                            radio = (size_point/5) # el 5 es constante

                            sam_circle_patches= {}
                            sam_circle_patches_pos = {}
                            sam_text_nodos = {}
                            for xx, yy, label, c in zip(DF_3.PCA1, DF_3.PCA2, id_samples, col):

                                if markers_point == '⏺':
                                    CC = ax.add_patch(Circle((xx, yy), radio, angle= 0, facecolor= c, zorder=2, alpha = 0.75, edgecolor = 'white', linewidth = 0.75,
                                                            url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[label]))
                                    sam_circle_patches[label] = CC

                                if markers_point == '⏹':
                                    CC = ax.add_patch(Rectangle((xx- (radio/2), yy - (radio/2)), radio, radio, angle = 0, facecolor=c, zorder=2, alpha = 0.75, edgecolor = 'white', linewidth=0.75,
                                                          url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[label]))
                                    sam_circle_patches[label] = CC


                                AA = ax.text(xx + (radio/2), yy, ' '+label.split('_')[-1], ha = 'left', va = 'center', color = 'black', fontsize = 8, alpha = 1, zorder = 2, family = '3ds Light')
                                sam_text_nodos[label] = AA
                                #ax.scatter(xx, yy, c= c, s = 150, marker = 'o', alpha = 0.8,  edgecolors='white', linewidths = 0.75, zorder=2)

                                sam_circle_patches_pos[label] = (xx, yy)


                            plt.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, zorder = 0)
                            plt.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, zorder = 0)

                            plt.xlabel('PCA 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)', size=LETRA_2pca, weight="bold")
                            plt.ylabel('PCA 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)', size=LETRA_2pca, weight="bold")


                            var_texto = {}
                            lines_plot = {}
                            
                            for cl in DF_3.clase.drop_duplicates():
                                points = []
                                a = DF_3[DF_3.clase == cl].PCA1.values
                                b = DF_3[DF_3.clase == cl].PCA2.values
                                #ax.scatter(np.mean(a), np.mean(b), zorder=0, c= 'white', edgecolors = dict_variable_element_colors[variable][cl], s = 3000, marker = 'o', alpha = 0.1, lw = 5)
                                GG = ax.text(np.mean(a), np.mean(b), cl, ha = 'center', va = 'center', color = 'black', fontsize = 15, alpha = 0.1, weight = 'bold', zorder = 2)
                                var_texto[cl] = GG

                                JJ = []
                                for A, B in zip(a, b):
                                    HH = plt.plot([np.mean(a), A], [np.mean(b), B],
                                             linestyle='--', markeredgewidth=0, zorder=0, markersize=0, color=category[cl], linewidth=1, alpha = 0.3)
                                    JJ.append(HH)

                                lines_plot[cl] = JJ


                            ax.text(ax.get_xlim()[0], ax.get_ylim()[1] + (ax.get_ylim()[1]*0.01), 'PERMANOVA: p='+str(Permanova)+'\n',
                                     fontsize=9, ha='left', va = 'center')


                            proxies = []
                            labels = []
                            for cat in category:
                                proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                                         c=category[cat], marker='o', alpha = 0.75,
                                                         markersize = 7,
                                                         markeredgecolor=category[cat], linewidth = 0)
                                proxies.append(proxy)
                                labels.append(cat)


                            ax.legend(proxies, list(category.keys()), title = ' '.join(["$\\bf{"+i+"}$" for i in VariablE.split(' ')]), numpoints=1, loc=2,borderpad = 0.2,
                                      handletextpad=-0.2,prop={'size':9},
                                              bbox_to_anchor=(1.005, 1.015))


                            ax.tick_params(bottom=True, right=False, top=False, left=True, width = 2, length=4, color='gainsboro')
                            #plt.gca().tick_params(which='major', width = 2, length=4, color='red')
                            for lados in ['left', 'right', 'top', 'bottom']:
                                plt.gca().spines[lados].set_linewidth(2)
                                plt.gca().spines[lados].set_color('gainsboro')


                            ax0 = fig.add_axes([1.3, 0, 0.3, 0.75]) 

                            ax0.set_xlim([dvw, -0.2])
                            ax0.set_ylim([0, ivw])

                            ax0.set_yticks(iv_ticks)

                            ax0.yaxis.set_ticks_position('right')


                            for line in ax0.get_yticklines():
                                line.set_visible(False)

                            colors_used = _remove_dups(color_list)
                            color_to_lines = {}
                            for color in colors_used:
                                color_to_lines[color] = []
                            for (xline, yline, color) in zip(xlines, ylines, color_list):
                                color_to_lines[color].append(list(zip(xline, yline)))

                            colors_to_collections = {}
                            for color in colors_used:
                                coll = matplotlib.collections.LineCollection(color_to_lines[color], colors=(color,))
                                colors_to_collections[color] = coll

                            for color in colors_used:
                                if color != above_threshold_color:
                                    ax0.add_collection(colors_to_collections[color])

                            if above_threshold_color in colors_to_collections:
                                    ax0.add_collection(colors_to_collections[above_threshold_color])

                            for e, i in enumerate(ax0.collections):
                                i.set_color(Set1[e])

                            ax0.text(0, ivw + (ivw*0.05), 'Clustering', fontsize=10, ha='left', va = 'center', weight = 'bold')

                            scatter_ax0 = {}
                            texto_ax0 = {}
                            for i, j in zip(ax0.get_yticks(), ivl):
                                OO = ax0.text(-0.13, i, j, ha = 'left', va = 'center', color = 'black', fontsize = 9, alpha = 1, zorder = 2)
                                texto_ax0[j] = OO
                                LL = ax0.scatter(-0.07, i, c= category[sam_clase[j]], s = 50, marker = 'o', alpha = 1,  edgecolors='none', linewidths = 0.75, zorder=2)
                                scatter_ax0[j] = LL

                            ax0.axis('off')

                            plt.close()


                            def set_locus():

                                if markers_point == '⏺':
                                    [sam_circle_patches[i].set_radius(size_point/5) for i in sam_circle_patches]

                                    proxies = []
                                    labels = []
                                    for cat in category:
                                        proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                                                 c=category[cat], marker='o', alpha = 0.75,
                                                                 markersize = 7,
                                                                 markeredgecolor=category[cat], linewidth = 0)
                                        proxies.append(proxy)
                                        labels.append(cat)


                                    ax.legend(proxies, list(category.keys()), title = ' '.join(["$\\bf{"+i+"}$" for i in VariablE.split(' ')]), numpoints=1, loc=2,borderpad = 0.2,
                                              handletextpad=-0.2,prop={'size':9, 'family':family_axis},
                                                      bbox_to_anchor=(1.005, 1.015))

                                if markers_point == '⏹':
                                    # cambia el tamano de los markers, pero al cambiar el tamano se camvia la posocion
                                    [sam_circle_patches[i].set_width((size_point/5)*1.7) for i in sam_circle_patches]
                                    [sam_circle_patches[i].set_height((size_point/5)*1.7) for i in sam_circle_patches]
                                    # con esto se reajustan las posiciones de los nuevos markers
                                    [sam_circle_patches[i].set_x(sam_circle_patches_pos[i][0] - ((size_point/5)/2)*1.7) for i in sam_circle_patches]
                                    [sam_circle_patches[i].set_y(sam_circle_patches_pos[i][1] - ((size_point/5)/2)*1.7) for i in sam_circle_patches]


                                    proxies = []
                                    labels = []
                                    for cat in category:
                                        proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                                                 c=category[cat], marker='s', alpha = 0.75,
                                                                 markersize = 7,
                                                                 markeredgecolor=category[cat], linewidth = 0)
                                        proxies.append(proxy)
                                        labels.append(cat)


                                    ax.legend(proxies, list(category.keys()), title = ' '.join(["$\\bf{"+i+"}$" for i in VariablE.split(' ')]), title_fontsize = 10, numpoints=1, loc=2,borderpad = 0.2,
                                              handletextpad=-0.2,prop={'size':9, 'family':family_axis},
                                                      bbox_to_anchor=(1.005, 1.015))


                                if ver_dendograma == 'False':
                                    #[scatter_ax0[i].set_alpha(0) for i in scatter_ax0]
                                    #[texto_ax0[i].set_alpha(0) for i in texto_ax0]
                                    #for e, i in enumerate(ax0.collections):
                                    #    i.set_linewidth(0)
                                    ax0.clear()
                                    ax0.axis('off')

                                if ver_dendograma == 'True':
                                    pass

                                if ver_num_muestra == 'False':
                                    [sam_text_nodos[i].set_alpha(0) for i in sam_text_nodos]
                                if ver_num_muestra == 'True':
                                    pass


                                [scatter_ax0[i].set_paths([marcas[{'⏺':'o', '⏹':'s'}[markers_point]]]) for i in scatter_ax0]

                                [var_texto[i].set_alpha(ver_var_name) for i in var_texto]

                                


                                [sam_circle_patches[i].set_alpha(alfa_point) for i in sam_circle_patches]

                                [j[0].set_alpha(alfa_linea) for i in lines_plot for j in lines_plot[i]]
                                
                                [texto_ax0[i].set_family(family_den) for i in texto_ax0]
                                
                                for labejex in ax.xaxis.get_ticklabels():
                                    labejex.set_fontsize(ticklabels_X)
                                for labejey in ax.yaxis.get_ticklabels():
                                    labejey.set_fontsize(ticklabels_Y)
                                for tickx in ax.xaxis.get_major_ticks():
                                    tickx.label.set_fontfamily([family_axis]) 
                                for ticky in ax.yaxis.get_major_ticks():
                                    ticky.label.set_fontfamily([family_axis])
                                ax.xaxis.get_label().set_fontsize(label_X)
                                ax.yaxis.get_label().set_fontsize(label_Y)
                                ax.xaxis.get_label().set_fontfamily([family_axis])
                                ax.yaxis.get_label().set_fontfamily([family_axis])
                                
                                if cambiar_sample == 'Code1':
                                    pass
                                if cambiar_sample == 'Code2':
                                    [texto_ax0[i].set_text(name_code[texto_ax0[i].get_text()]) for i in texto_ax0]
                                if cambiar_sample == 'Code3':
                                    [texto_ax0[i].set_text(name_code2[texto_ax0[i].get_text()]) for i in texto_ax0]


                                display(fig)

                        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                        if select_PCA == '3D':

                            plot_alto_pca3 = 9
                            LETRA_3pca = 9

                            mpl.rcParams.update(mpl.rcParamsDefault)
                            plt.rcParams['grid.color'] = "whitesmoke"


                            fig = plt.figure(figsize=(plot_alto_pca3, plot_alto_pca3-2))

                            ax = fig.add_subplot(111, projection = '3d', facecolor = 'white')

                            ax.grid(True)
                            #ax.xaxis.pane.fill = False
                            #ax.yaxis.pane.fill = False
                            #ax.zaxis.pane.fill = False
                            
                            sam_text_nodos = {}
                            sam_circle_patches = {}
                            for x, y, z, label, c  in zip(DF_3.PCA1, DF_3.PCA2, DF_3.PCA3, id_samples, col):
                                FF = ax.scatter(x, y, z, c= c, s = 100, marker = 'o', alpha = 0.7, edgecolors='white', linewidths = 1, zorder=1)
                                sam_circle_patches[label] = FF
                                
                                
                                AA = ax.text(x, y, z, '   '+label.split('_')[-1], ha = 'left', va = 'center', color = 'black', fontsize = 9, alpha = 1, zorder = 2, family = '3ds Light')
                                sam_text_nodos[label] = AA
                                

                            ax.text(ax.get_xlim()[1] * 0.1, ax.get_ylim()[1], ax.get_zlim()[1], 'PERMANOVA: p='+str(Permanova)+'\n\n\n\n',
                                    fontsize=9, ha='left', va = 'center')

                            proxies = []
                            labels = []
                            for cat in category:
                                proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                                         c=category[cat], marker={'⏺':'o', '⏹':'s'}[markers_point], alpha = 0.7,
                                                         markersize = 7,
                                                         markeredgecolor=category[cat], linewidth = 0)
                                proxies.append(proxy)
                                labels.append(cat)


                            ax.legend(proxies, list(category.keys()), title = ' '.join(["$\\bf{"+i+"}$" for i in variable.split(' ')]), numpoints=1, loc=2,borderpad = 0.2,
                                      handletextpad=-0.2,prop={'size':9, 'family':family_axis},
                                              bbox_to_anchor=(1.015, 0.93))
                            
                            var_texto = {}
                            lines_plot = {}
                            for cl in DF_3.clase.drop_duplicates():
                                points = []
                                a = DF_3[DF_3.clase == cl].PCA1.values
                                b = DF_3[DF_3.clase == cl].PCA2.values
                                c = DF_3[DF_3.clase == cl].PCA3.values
                                #ax.scatter(np.mean(a), np.mean(b), np.mean(c), zorder=0, c= 'black', s = 5, marker = 'o', alpha = 1, lw = 0)
                                GG = ax.text(np.mean(a), np.mean(b), np.mean(c), cl, ha = 'center', va = 'center', color = 'black', fontsize = 15, alpha = 0.1, weight = 'bold', zorder = 2)
                                var_texto[cl] = GG
                                
                                JJ = []
                                for A, B, C in zip(a, b, c):
                                    HH = plt.plot([np.mean(a), A], [np.mean(b), B], [np.mean(c), C],
                                             linestyle='--', markeredgewidth=0, zorder=0, markersize=0, color=category[cl], linewidth=0.5, alpha = 0.3)
                                    JJ.append(HH)
                                lines_plot[cl] = JJ

                            ax.set_xlabel('PCA 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)', size=LETRA_3pca, weight="bold")
                            ax.set_ylabel('PCA 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)', size=LETRA_3pca, weight="bold")
                            ax.set_zlabel('PCA 3 ('+str(round(pca3.explained_variance_ratio_[2]*100,2))+'%)', size=LETRA_3pca, weight="bold")

                            plt.xticks(size=9)
                            plt.yticks(size=9)

                            ax.zaxis.set_tick_params(labelsize=9)

                            fccolor3D = 'white'
                            ax.xaxis.pane.set_facecolor(fccolor3D)
                            ax.yaxis.pane.set_facecolor(fccolor3D)
                            ax.zaxis.pane.set_facecolor(fccolor3D)
                            ax.xaxis.pane.set_edgecolor('gainsboro')
                            ax.yaxis.pane.set_edgecolor('gainsboro')
                            ax.zaxis.pane.set_edgecolor('gainsboro')
                            ax.xaxis.pane.set_alpha(1)
                            ax.yaxis.pane.set_alpha(1)
                            ax.zaxis.pane.set_alpha(1)
                            ax.xaxis.pane.set_linewidth(2)
                            ax.yaxis.pane.set_linewidth(2)
                            ax.zaxis.pane.set_linewidth(2)

                            #---------------------------------------------

                            ax0 = fig.add_axes([1.18, 0.1, 0.2, 0.65]) 

                            ax0.set_xlim([dvw, -0.2])
                            ax0.set_ylim([0, ivw])

                            ax0.set_yticks(iv_ticks)

                            ax0.yaxis.set_ticks_position('right')


                            for line in ax0.get_yticklines():
                                line.set_visible(False)

                            colors_used = _remove_dups(color_list)
                            color_to_lines = {}
                            for color in colors_used:
                                color_to_lines[color] = []
                            for (xline, yline, color) in zip(xlines, ylines, color_list):
                                color_to_lines[color].append(list(zip(xline, yline)))

                            colors_to_collections = {}
                            for color in colors_used:
                                coll = matplotlib.collections.LineCollection(color_to_lines[color], colors=(color,))
                                colors_to_collections[color] = coll

                            for color in colors_used:
                                if color != above_threshold_color:
                                    ax0.add_collection(colors_to_collections[color])

                            if above_threshold_color in colors_to_collections:
                                    ax0.add_collection(colors_to_collections[above_threshold_color])

                            for e, i in enumerate(ax0.collections):
                                i.set_color(Set1[e])

                            ax0.text(0, ivw + (ivw*0.05), 'Clustering', fontsize=10, ha='left', va = 'center', weight = 'bold')
                            
                            scatter_ax0 = {}
                            texto_ax0 = {}
                            for i, j in zip(ax0.get_yticks(), ivl):
                                OO = ax0.text(-0.13, i, j, ha = 'left', va = 'center', color = 'black', fontsize = 9, alpha = 1, zorder = 2)
                                texto_ax0[j] = OO
                                LL = ax0.scatter(-0.07, i, c= category[sam_clase[j]], s = 50, marker = 'o', alpha = 1,  edgecolors='none', linewidths = 0.75, zorder=2)
                                scatter_ax0[j] = LL
                            
                            ax0.axis('off')

                            plt.close()

                            def set_locus():
                                

                                [j[0].set_alpha(alfa_linea) for i in lines_plot for j in lines_plot[i]]
                                
                                [sam_circle_patches[i].set_paths([marcas[{'⏺':'o', '⏹':'s'}[markers_point]]]) for i in sam_circle_patches]
                                
                                [scatter_ax0[i].set_paths([marcas[{'⏺':'o', '⏹':'s'}[markers_point]]]) for i in scatter_ax0]
                                
                                [texto_ax0[i].set_family(family_den) for i in texto_ax0]
                                
                                [sam_circle_patches[i].set_sizes([(size_point*2)*100]) for i in sam_circle_patches]
                                
                                [sam_circle_patches[i].set_alpha(alfa_point) for i in sam_circle_patches]
                                
                                [var_texto[i].set_alpha(ver_var_name) for i in var_texto]
                                
                                if cambiar_sample == 'Code1':
                                    pass
                                if cambiar_sample == 'Code2':
                                    [texto_ax0[i].set_text(name_code[texto_ax0[i].get_text()]) for i in texto_ax0]
                                if cambiar_sample == 'Code3':
                                    [texto_ax0[i].set_text(name_code2[texto_ax0[i].get_text()]) for i in texto_ax0]
                                
                                if ver_dendograma == 'False':
                                    ax0.clear()
                                    ax0.axis('off')

                                if ver_dendograma == 'True':
                                    pass
                                
                                for labejex in ax.xaxis.get_ticklabels():
                                    labejex.set_fontsize(ticklabels_X)
                                for labejey in ax.yaxis.get_ticklabels():
                                    labejey.set_fontsize(ticklabels_Y)
                                for labejey in ax.zaxis.get_ticklabels():
                                    labejey.set_fontsize(ticklabels_Z)
                                    
                                    
                                for tickx in ax.xaxis.get_major_ticks():
                                    tickx.label.set_fontfamily([family_axis]) 
                                for ticky in ax.yaxis.get_major_ticks():
                                    ticky.label.set_fontfamily([family_axis])
                                for ticky in ax.zaxis.get_major_ticks():
                                    ticky.label.set_fontfamily([family_axis])    
                                    
                                ax.xaxis.get_label().set_fontsize(label_X)
                                ax.yaxis.get_label().set_fontsize(label_Y)
                                ax.zaxis.get_label().set_fontsize(label_Z)
                                
                                ax.xaxis.get_label().set_fontfamily([family_axis])
                                ax.yaxis.get_label().set_fontfamily([family_axis])
                                ax.zaxis.get_label().set_fontfamily([family_axis])
                                
                                if ver_num_muestra == 'False':
                                    [sam_text_nodos[i].set_alpha(0) for i in sam_text_nodos]
                                if ver_num_muestra == 'True':
                                    pass
                                

                                display(fig)
                        
                        
                        # --- save
                        ancho_plot_save = str(71)
                        png1 = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
                        png1.style.button_color = 'gold'
                        output1 = widgets.Output()
                        def button_clicked1(b):
                            with output1:
                                clear_output(True)
                                if filename_plot.value is '':
                                    nombre_grafico = 'Beta_diversity_no_phyl_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                                if filename_plot.value is not '':
                                    nombre_grafico = filename_plot.value
                                    
                                fig.savefig('plots_asv/beta_diversity/'+nombre_grafico+'.png', dpi = 900, bbox_inches= 'tight')
                        
                        png1.on_click(button_clicked1)
                        #----
                        jpeg1 = widgets.Button(description="JPEG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
                        jpeg1.style.button_color = 'gold'
                        output2 = widgets.Output()
                        def button_clicked2(b):
                            with output2:
                                clear_output(True)
                                if filename_plot.value is '':
                                    nombre_grafico = 'Beta_diversity_no_phyl_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                                if filename_plot.value is not '':
                                    nombre_grafico = filename_plot.value
                                    
                                fig.savefig('plots_asv/beta_diversity/'+nombre_grafico+'.jpeg', dpi = 900, bbox_inches= 'tight')
                                
                                
                        jpeg1.on_click(button_clicked2)
                        #----
                        svg1 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
                        svg1.style.button_color = 'gold'
                        output3 = widgets.Output()
                        def button_clicked3(b):
                            with output3:
                                clear_output(True)
                                if filename_plot.value is '':
                                    nombre_grafico = 'Beta_diversity_no_phyl_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                                if filename_plot.value is not '':
                                    nombre_grafico = filename_plot.value
                                    
                                fig.savefig('plots_asv/beta_diversity/'+nombre_grafico+'.svg', dpi = 900, bbox_inches= 'tight')
                                
                                
                        svg1.on_click(button_clicked3)
                        #----
                        pdf1 = widgets.Button(description="PDF", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
                        pdf1.style.button_color = 'gold'
                        output4 = widgets.Output()
                        def button_clicked4(b):
                            with output4:
                                clear_output(True)
                                if filename_plot.value is '':
                                    nombre_grafico = 'Beta_diversity_no_phyl_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                                if filename_plot.value is not '':
                                    nombre_grafico = filename_plot.value
                                    
                                fig.savefig('plots_asv/beta_diversity/'+nombre_grafico+'.pdf', dpi = 900, bbox_inches= 'tight')
                                
                                
                        pdf1.on_click(button_clicked4)
                        ### ---


                        params = widgets.Button(description="SAVE", icon = 'fa-save', layout=Layout(width=ancho_plot_save+'px'))
                        params.style.button_color = 'cyan'
                        output5 = widgets.Output()
                        def button_clicked5(b):
                            with output5:
                                clear_output(True)

                                if filename.value is '':
                                    nombre_archivo = 'Beta_diversity_no_phyl_chart_params_ASVs.txt'
                                if filename.value is not '':
                                    nombre_archivo = filename.value

                                with open('plots_asv/beta_diversity/'+nombre_archivo, 'w') as fq:
                                    fq.write('#Saved parameters\n')
                                    fq.write('Clustering:ASVs\n' )
                                    """
                                    diccionario para guardar en un archivo los parametros usados para volver a reproducirlos posteriormente
                                    """
                                    parametros = {'Kit:':tipo_kits_1,
                                                  'Variable:':VariablE,
                                                  'Chart type:':select_PCA,
                                                  'Marker:':{'⏺':'circle', '⏹':'square'}[markers_point],
                                                  'Sample number:':ver_num_muestra,
                                                  'Marker size:':size_point,
                                                  'Marked alpha:':alfa_point,
                                                  'Line alpha:':alfa_linea,
                                                  'Text alpha:':ver_var_name,
                                                  
                                                  'Dendrogram font:':family_den,
                                                  'Dendrogram:':ver_dendograma,
                                                  'Rename sample:':cambiar_sample,
                                                  'X label:':label_X,
                                                  'X tick label:':ticklabels_X,
                                                  'Y label:':label_Y,
                                                  'Y tick labels:':ticklabels_Y,
                                                  'Z label:':label_Y,
                                                  'Z tick labels:':ticklabels_Y,
                                                  'Axis font:':family_axis,
                                                  '\n###\nWarning:':'The variables take the same colors assigned in the taxonomy analysis.'}
                                    for w in parametros:
                                        fq.write(w+str(parametros[w])+'\n')

                        params.on_click(button_clicked5)
                        #------
                        
                        items_save_1 = VBox([HBox([widgets.HTML('<font color = grey> <b style="font-size:0.7vw">SAVE CHART: </b>'), blanca, filename_plot]),
                                             HBox([blanca, widgets.Label('Formats:'), png1, jpeg1, svg1, pdf1]),
                                     HBox([blanca, widgets.Label('Chart parameters:'), filename, params])])
                        items_save_1_box = Box(children=[items_save_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_save_1.children) * 34))+'px'))
                        
                        display(VBox([HBox([RAMPA_ELEGIDA, items_save_1_box]), output5]))
                        
                        set_locus()
                        

                    OUT_BetA2 = widgets.interactive_output(BetA2, {'markers_point':markers_point, 'size_point':size_point, 'ver_num_muestra':ver_num_muestra,
                                                     'ver_dendograma':ver_dendograma, 'ver_var_name':ver_var_name, 
                                                     'alfa_point':alfa_point, 'alfa_linea':alfa_linea,
                                                                  'label_X':label_X, 'ticklabels_X':ticklabels_X, 'label_Y':label_Y, 'ticklabels_Y':ticklabels_Y, 'family_axis':family_axis,
                                                                   'label_Z':label_Z, 'ticklabels_Z':ticklabels_Z, 'family_den':family_den, 'cambiar_sample':cambiar_sample,
                                                                  'SalidA0':SalidA0})
                    
                    
                    display(OUT_BetA2)


                if select_PCA == '3Di':
                    import datetime
                    category = dict_variable_element_colors[VariablE]
                    
                    blanca0 = widgets.Button(layout=Layout(width='3px', height='25px'), disabled=True)
                    blanca0.style.button_color = 'white'

                    # --- save
                    html1 = widgets.Button(description="HTML", icon = 'fa-bar-chart', layout=Layout(width='71px'))
                    html1.style.button_color = 'gold'
                    output11 = widgets.Output()
                    def button_clicked11(b):
                        with output11:
                            clear_output(True)
                            plotly.offline.plot(FigurA, filename = 'plots_asv/beta_diversity/3Di_'+tipo_kits_1+'_'+re.sub('[/]|[(]|[)]|[*]|[|]', '_', VariablE)+'_ASVs'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.html', auto_open=False)

                    html1.on_click(button_clicked11)

                    ver_linea = widgets.ToggleButtons(options=['True', 'False'], value = 'True')
                    ver_linea.style.button_width = '55px'
                    ver_texto = widgets.ToggleButtons(options=['True', 'False'], value = 'True')
                    ver_texto.style.button_width = '55px'
                    ver_estilo = widgets.ToggleButtons(options=['White', 'Black'], value = 'White')
                    ver_estilo.style.button_width = '55px'


                    desc_ver_linea = HBox([widgets.Label('Show line:'), ver_linea, blanca0, widgets.Label('Sample number:'), ver_texto, blanca0, widgets.Label('Style:'), ver_estilo, blanca0, widgets.Label('Save 3Di:'), html1])
                    plot_3di = HBox([widgets.HTML('<font color = black> <b style="font-size:0.6vw">3Di Settings: </b>'), blanca0, desc_ver_linea])
                    set_box = Box(children=[plot_3di], layout=Layout(border='1px solid black', width='898px', height='37px'))



                    Lines_poS = {}
                    for cl in DF_3.clase.drop_duplicates():
                        a = DF_3[DF_3.clase == cl].PCA1.values
                        b = DF_3[DF_3.clase == cl].PCA2.values
                        c = DF_3[DF_3.clase == cl].PCA3.values

                        lines_positions_data = {}
                        n = 0
                        for A, B, C in zip(a, b, c):
                            n += 1        
                            lines_positions_data['data'+str(n)] = {'x':np.array([np.mean(a), A]), 'y':np.array([np.mean(b), B]), 'z':np.array([np.mean(c), C])}
                        Lines_poS[cl] = lines_positions_data


                    DF3D = DF_3
                    DF3D['Size'] = 15
                    DF3D[VariablE] = DF3D.clase
                    DF3D['Sample'] = DF3D.index
                    DF3D = DF3D.reset_index(drop = True)

                    colorr = list(category.values())

                    labels = {'PCA1': 'PC 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)',
                              'PCA2': 'PC 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)',
                              'PCA3': 'PC 3 ('+str(round(pca3.explained_variance_ratio_[2]*100,2))+'%)'}

                    fig = go.Figure()

                    scatter_clases = 0
                    for i in DF3D.clase.unique():
                        x = DF3D[DF3D.clase == i].PCA1.values
                        y = DF3D[DF3D.clase == i].PCA2.values
                        z = DF3D[DF3D.clase == i].PCA3.values
                        lab = DF3D[DF3D.clase == i].Sample.tolist()
                        fig.add_trace(go.Scatter3d(x=x, y=y, z=z,
                                                           name = i,
                                                           projection=None,
                                                           text= 'none',
                                                           hovertext = [VariablE+'='+i+'<br>'+labels['PCA1']+'='+str(round(a, 5))+'<br>'+labels['PCA2']+'='+str(round(b, 5))+'<br>'+labels['PCA3']+'='+str(round(c, 5))+'<br>'+k for k, a, b, c in zip(lab, x, y, z)],
                                                           hoverinfo='text',
                                                           marker_color = category[i],
                                                           opacity=1,
                                                           showlegend=True,
                                                           marker=dict(size=6, opacity=1, line=dict(color=category[i], width=200)),
                                                           textfont=dict(family="3ds Light", size=10, color=category[i]),
                                                           mode='markers'))
                        scatter_clases += 1


                    scatter_lines = scatter_clases
                    for h in Lines_poS:
                        for i in Lines_poS[h]:
                            fig.add_trace(go.Scatter3d(x=Lines_poS[h][i]['x'], y=Lines_poS[h][i]['y'], z=Lines_poS[h][i]['z'], hoverinfo='none', showlegend=False,
                                                        opacity=0.5, mode='lines', name='lines', line = dict(color=category[h], width=5, dash='solid')))
                            scatter_lines += 1


                    scater_con_text = scatter_lines
                    for i in DF3D.clase.unique():
                        x = DF3D[DF3D.clase == i].PCA1.values
                        y = DF3D[DF3D.clase == i].PCA2.values
                        z = DF3D[DF3D.clase == i].PCA3.values
                        lab = DF3D[DF3D.clase == i].Sample.tolist()
                        fig.add_trace(go.Scatter3d(x=x, y=y, z=z,
                                                           name = i,
                                                           text = [i.split('_')[-1] for i in lab],
                                                           projection=None,
                                                           hovertext = [VariablE+'='+i+'<br>'+labels['PCA1']+'='+str(round(a, 5))+'<br>'+labels['PCA2']+'='+str(round(b, 5))+'<br>'+labels['PCA3']+'='+str(round(c, 5))+'<br>'+k for k, a, b, c in zip(lab, x, y, z)],
                                                           hoverinfo='text',
                                                           marker_color = category[i],
                                                           opacity=1,
                                                           showlegend=False,
                                                           marker=dict(size=6, opacity=1, line=dict(color=category[i], width=200)),
                                                           textfont=dict(family="3ds Light", size=13),
                                                           mode='markers+text'))
                        scater_con_text += 1


                    fig.update_layout(legend_title_text=VariablE, autosize=True,width=900, height=650,margin=dict(l=0, r=5, b=0, t=0, pad=0),
                                      scene = dict(xaxis_title=labels['PCA1'], yaxis_title=labels['PCA2'], zaxis_title=labels['PCA3']),
                                         font_size=13, template="plotly_white", legend=dict(yanchor="top", y=0.9, xanchor="left", x=1))

                    FigurA = go.FigureWidget(fig)
                    display(VBox([set_box, FigurA]))

                    def res(f):
                        with FigurA.batch_update():
                            if ver_linea.value == 'True':
                                for i in FigurA.data[scatter_clases:scatter_lines]:
                                    i.opacity = 0.5
                            if ver_linea.value == 'False':
                                for i in FigurA.data[scatter_clases:scatter_lines]:
                                    i.opacity = 0

                            if ver_texto.value == 'True':
                                for i in FigurA.data[scatter_lines:]:
                                    i.mode = 'markers+text'
                            if ver_texto.value == 'False':
                                for i in FigurA.data[scatter_lines:]:
                                    i.mode = 'markers'
                                    
                            if ver_estilo.value == 'White':
                                FigurA.layout.template = "plotly_white"
                                
                            if ver_estilo.value == 'Black':
                                FigurA.layout.template = "plotly_dark"


                    ver_linea.observe(res, names="value")
                    ver_texto.observe(res, names="value")
                    ver_estilo.observe(res, names="value")
                
        OUT_BetA = widgets.interactive_output(BetA, {'tipo_kits_1':tipo_kits_1, 'VariablE':VariablE, 'select_PCA':select_PCA,})
        
        items_data = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'DATA']
        items_axis = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'AXIS']
        items_save = [widgets.HTML('<font color = white> <b style="font-size:0.7vw">'+i+'</b>') for i in 'S']
        
        
        
        
        items_plot = VBox([HBox([blanca, widgets.Label('Marker size:'),
                                 markers_point,
                                 widgets.Label('Sample number:'), ver_num_muestra]),
                                 size_point,
                                 alfa_point,
                                 alfa_linea,
                                 ver_var_name,
                                 HBox([widgets.Label('Dendrogram font:'), family_den]),
                                 HBox([widgets.Label('Dendrogram:'), ver_dendograma]),
                                 HBox([widgets.Label('Rename sample:'), cambiar_sample])])
        
        items_plot_box = Box(children=[items_plot], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_plot.children) * 31))+'px'))
        
        items_axis_1 = VBox([label_X,
                                 ticklabels_X,
                                 label_Y,
                                 ticklabels_Y,
                                 label_Z,
                                 ticklabels_Z,
                                 HBox([blanca, widgets.Label('Axis font:'), family_axis])])
        
        items_axis_1_box = Box(children=[items_axis_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_axis_1.children) * 31))+'px'))
        
        
        display(HBox([VBox([VBox([widgets.HTML('<font color = #1976d2> <b style="font-size:0.8vw">KITS</b>'), tipo_kit_1_box]),
                            VBox([widgets.HTML('<font color = gray> <b style="font-size:0.8vw">VARIABLES</b>'), tipo_variable_box]),
                            VBox([widgets.HTML('<font color = gray> <b style="font-size:0.6vw">Chart type:</b>'), boton_select_PCA])]), blanca,
                      VBox([widgets.HTML('<font color = grey> <i class="fa fa-cog fa-2x fa-fw"></i> <b style="font-size:0.8vw">PLOT SETTINGS (2D & 3D)</b>'),
                            HBox([Box(children =[VBox(items_data)]), items_plot_box]),
                             
                            HBox([Box(children =[VBox(items_axis)]), items_axis_1_box]),
                            
                           
                           ]),
                      blanca, OUT_BetA]))
        
BETA_button.on_click(button_clicked)





BETA_DIVERSITY_NO_PHYL = VBox([HBox([BETA_button, estatico_no_fil_box]), BETA_output])




















# # BETA DIVERSITY PHYLOGENETIC




from io import StringIO
from skbio.tree import TreeNode





# dos niveles
import datetime

progreso1_fil = widgets.HTML('<font color = black> <i class="fa fa-spinner fa-2x fa-fw"></i> </font>')
progreso2_fil = widgets.HTML()
estatico_fil = widgets.HTML('<font color = limegreen> <b style="font-size:0.5vw">Processed samples : </b>')
estatico_fil = HBox([progreso1_fil, estatico_fil, progreso2_fil])
estatico_fil_box = Box(children=[estatico_fil], layout=Layout(border='1px solid limegreen', width='585px', height='32px'))

BETA_PHY_button = widgets.Button(description="PROCESS AND VISUALIZE", icon = 'fa-eye', layout=Layout(width='590px'))
BETA_PHY_button.style.button_color = 'gainsboro' #'deepskyblue'
BETA_PHY_button.style.font_weight = 'bold'
BETA_PHY_output = widgets.Output()

def button_clicked(b):
    import time
    from scipy.cluster.hierarchy import dendrogram, average
    
    
    NCBI_RDP_SILVA_SUMMARY = pd.read_csv('tablas/ASVs_NCBI_RDP_SILVA_SUMMARY.txt', sep = '\t')
    
    rarefaction_collapse_specie = pd.pivot_table(NCBI_RDP_SILVA_SUMMARY[['Species'] + name_sample], values = name_sample, index = ['Species'], aggfunc = sum).reset_index()
    
    with open('plots_asv/taxonomy/dict_variable_element_colors_ALL.json', 'r') as fp:
        dict_variable_element_colors = json.load(fp)

    asv_species_dict = {}
    asv_species_lista = []
    for i in NCBI_RDP_SILVA_SUMMARY.Species.drop_duplicates():
        w = NCBI_RDP_SILVA_SUMMARY[NCBI_RDP_SILVA_SUMMARY.Species == i]
        sumas = np.sum(w.iloc[:, 3:].values, axis = 1)
        #sumas = np.sum(w.values[0].tolist()[3:])
        w['sum'] = sumas
        w = w.sort_values(by =['sum'],ascending=False).reset_index(drop=True)
        ef = w['Entry'].tolist()[0]
        asv_species_dict[ef] = i
        asv_species_lista.append([ef, i])
    unicos = DataFrame(asv_species_lista, columns = ['Entry', 'Species'])
    
    MERGE_UNICOS_SPECIES = unicos.merge(rarefaction_collapse_specie, on = 'Species', how = 'left')
    
    arbol1 = open('clustering/ASVs_Phylogenetic_Diversity_tree.txt', 'r')
    arbol = arbol1.read()
    arbol1.close()
    arbol = ''.join(arbol.split('\n'))
    arbol = re.sub(';$', 'root;', arbol)
    
    # este arbol filogenetico lo obtuve desde clustal omega
    newick_tree = StringIO(arbol)
    
    tree = TreeNode.read(newick_tree)
    tree = tree.root_at_midpoint()
    
    def get_observed_nodes(tree, table, sample_id, columna, verbose=False):

        observed_otus = table[[sample_id, columna]][table[sample_id] > 0][columna].tolist()

        #observed_otus = [obs_id for obs_id in table.index if table[sample_id][obs_id] > 0]

        observed_nodes = set()
        # iterate over the observed OTUs
        for otu in observed_otus:
            t = tree.find(otu)
            observed_nodes.add(t)
            if verbose:
                print(t.name, t.length, end=' ')
            for internal_node in t.ancestors():
                if internal_node.length is None:
                    # we've hit the root
                    if verbose:
                        pass
                else:
                    if verbose and internal_node not in observed_nodes:
                        print(internal_node.length, end=' ')
                    observed_nodes.add(internal_node)
        return observed_nodes
    def unweighted_unifrac(tree, table, sample_id1, sample_id2, columna, verbose=False):
        observed_nodes1 = get_observed_nodes(tree, table, sample_id1, columna, verbose=verbose)
        observed_nodes2 = get_observed_nodes(tree, table, sample_id2, columna, verbose=verbose)
        observed_branch_length = sum(o.length for o in observed_nodes1 | observed_nodes2)
        shared_branch_length = sum(o.length for o in observed_nodes1 & observed_nodes2)
        unique_branch_length = observed_branch_length - shared_branch_length
        unweighted_unifrac = unique_branch_length / observed_branch_length
        return unweighted_unifrac
    def table_to_distances2(table, columna, pairwise_distance_fn):
        sample_ids = table.iloc[:, 1:].columns
        num_samples = len(sample_ids)
        data = zeros((num_samples, num_samples))
        for i, sample1_id in enumerate(sample_ids):
            for j, sample2_id in enumerate(sample_ids[:i]):
                data[i,j] = data[j,i] = pairwise_distance_fn(tree, table, sample1_id, sample2_id, columna, verbose=False)
        return data, sample_ids.tolist()

    
    with BETA_PHY_output:
        clear_output(True)
        
        progreso1_fil.value = '<font color = black> <i class="fa fa-spinner fa-pulse fa-2x fa-fw"></i> </font>'
        
        
        PCA_KITS = {}
        for k in KITS:
            if k == 'Both kits':
                NEW_rarefaction_collapse_specie = MERGE_UNICOS_SPECIES[['Entry'] + KITS['Both kits']]
            if k == 'DNPowerSoil':
                NEW_rarefaction_collapse_specie = MERGE_UNICOS_SPECIES[['Entry'] + KITS['DNPowerSoil']]
            if k == 'DNMicrobial':
                NEW_rarefaction_collapse_specie = MERGE_UNICOS_SPECIES[['Entry'] + KITS['DNMicrobial']]
            
            MatriX, matrixLabels = table_to_distances2(NEW_rarefaction_collapse_specie, 'Entry', unweighted_unifrac)
            
            
            SAMple = DataFrame(matrixLabels, columns = ['Name Sample']).merge(metadata, on = 'Name Sample', how = 'left')
            SAMple =SAMple.set_index('Name Sample')
            
            MatriX_unweighted = DataFrame(MatriX, columns = matrixLabels)
            MatriX_unweighted.insert(loc = 0, column='Name Sample', value=matrixLabels)
            MatriX_unweighted = MatriX_unweighted.merge(metadata, on = 'Name Sample', how = 'left')
            MatriX_unweighted = MatriX_unweighted.set_index('Name Sample')
            
            
            SS = StandardScaler()
            MatriX_unweighted[matrixLabels] = SS.fit_transform(MatriX_unweighted[matrixLabels].astype('float'))

            pca3 = PCA(n_components=3)
            pca_3 = pca3.fit_transform(MatriX_unweighted[matrixLabels].astype('float'))
            

            pca_kits_variables = {}

            for variable in variables:
                progreso2_fil.value = '<font color = red> <b style="font-size:0.5vw">'+k+' = '+variable+'</b>'
                
                if len(MatriX_unweighted[variable].unique()) == 1:
                    pass
                else:
                    DF_3 = pd.DataFrame({'PCA1':pca_3[:,0], 'PCA2':pca_3[:, 1], 'PCA3':pca_3[:, 2], 'clase':MatriX_unweighted[variable]})
                    Permanova = permanova(MatriX, SAMple, variable)['p-value']

                    bc_matrix_condensed_ASVs = squareform(MatriX, force='tovector', checks=False)
                    ZZZ = average(bc_matrix_condensed_ASVs)

                    dendro = dendrogram(ZZZ,
                                orientation='left',  no_plot = True,
                                labels=matrixLabels,
                                distance_sort='descending',
                                        leaf_font_size = 100,
                                show_leaf_counts=True)
                    #-----------------------------------------
                    ivl = dendro['ivl']
                    color_list = dendro['color_list']
                    Z = np.asarray(dendro['dcoord'], order='c')
                    mh = max(Z[:, 2])
                    above_threshold_color = 'b'
                    # Independent variable plot width
                    ivw = len(ivl) * 10
                    # Dependent variable plot height
                    dvw = mh + mh * 0.05
                    iv_ticks = np.arange(5, len(ivl) * 10 + 5, 10)
                    xlines = dendro['dcoord']
                    ylines = dendro['icoord']

                    category = dict_variable_element_colors[variable]
                    col = DF_3['clase'].map(category)

                pca_kits_variables[variable] = {'df_3':DF_3, 'permanova':Permanova, 'ivl':ivl, 'ivw':ivw, 'dvw':dvw, 'category':category,
                                                'iv_ticks':iv_ticks, 'xlines':xlines, 'ylines':ylines, 'color_list':color_list, 'above_threshold_color':above_threshold_color}
                time.sleep(0.1)
                
                
            PCA_KITS[k] = pca_kits_variables
            
        del DF_3, Permanova, ivl, ivw, dvw, category, iv_ticks, xlines, ylines, color_list, above_threshold_color, matrixLabels
        
        clear_output(True)
        progreso1_fil.value = '<font color = black> <i class="fa fa-spinner fa-2x fa-fw"></i> </font>'
        progreso2_fil.value = '<font color = red> <b style="font-size:0.5vw">'+k+' = '+variable+'</b>'
        
        
        mpl.rcParams.update(mpl.rcParamsDefault)

        fig = plt.figure(figsize=(6, 4))

        ax = fig.add_axes([0, 0, 1, 1])

        plt.gca().tick_params(which='major', width = 2, length=4, color='gainsboro')
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_color('gainsboro')
        plt.gca().spines['bottom'].set_color('gainsboro')
        plt.gca().spines['right'].set_color(None)
        plt.gca().spines['top'].set_color(None)

        plt.gca().set_xlabel('Contigs', fontsize=12, fontname='Open Sans', weight = 'bold')
        plt.gca().set_ylabel('Species', fontsize=12, fontname='Open Sans', weight = 'bold')
        plt.close()
        
        
        blanca = widgets.Button(layout=Layout(width='10px', height='25px'), disabled=True)
        blanca.style.button_color = 'white'
        blanca2 = widgets.Button(layout=Layout(width='1px', height='25px'), disabled=True)
        blanca2.style.button_color = 'white'
        
        tipo_kits_1 = widgets.ToggleButtons(options= ['Both kits'] + metadata[VARIABLE_KIT].unique().tolist(), value = 'Both kits', button_style = 'primary')
        tipo_kits_1.style.button_width = '170px'
        tipo_kits_1.style.font_weight = 'bold'
        tipo_kit_1_box = Box(children=[VBox([tipo_kits_1])], layout=Layout(border='1px solid #1976d2', width='180px', height='95px'))

        VariablE = widgets.ToggleButtons(options= variables, value = variables[0], button_style = 'warning')
        VariablE.style.button_width = '170px'
        tipo_variable_box = Box(children=[VBox([VariablE])], layout=Layout(border='1px solid #ff9800', width='180px', height='210px'))
        

        select_PCA = widgets.ToggleButtons(options=['2D', '3D', '3Di'], value = '2D')
        select_PCA.style.button_width = '50px'
        boton_select_PCA = Box(children=[select_PCA], layout= Layout(border='1px solid pink', width='170px', height='35px'))
        
        markers_point = widgets.ToggleButtons(options=['⏺', '⏹'], value = '⏺')
        markers_point.style.button_width = '31px'
        
        size_point = widgets.SelectionSlider(options=np.round(np.arange(0, 1.505, 0.05), 2), value=0.5,disabled=False,
                                              description = 'Marker size:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        
        alfa_point = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.7,disabled=False,
                                              description = 'Marker alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        alfa_linea = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.3,disabled=False,
                                              description = 'Line alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px')) 
        
        ver_num_muestra = widgets.ToggleButtons(options=['True', 'False'], value = 'False')
        ver_num_muestra.style.button_width = '55px'
        
        
        ver_dendograma = widgets.ToggleButtons(options=['True', 'False'], value = 'True')
        ver_dendograma.style.button_width = '55px'
        
        ver_var_name = widgets.SelectionSlider(options=np.round(np.linspace(0, 1, 11), 2),value=0.1,disabled=False,
                                              description = 'Text alpha:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                           layout=Layout(width='400px', height='25px'))
        
        
        family = sorted(['Liberation Serif','Microsoft Sans Serif','Open Sans','Times New Roman','3ds Light','Calibri','Comic Sans MS',
                  'Arial','Courier New','Microsoft Yi Baiti','Lucida Console'])
        family_den = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                   layout = Layout(width='290px', height='25px'))

        ### ejes
        label_X = widgets.SelectionSlider(options=range(5, 31),value=12,disabled=False,
                                              description = 'X label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_X = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'X tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))

        label_Y = widgets.SelectionSlider(options=range(5, 31),value=12,disabled=False,
                                              description = 'Y label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        ticklabels_Y = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'Y tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))

        family_axis = widgets.Dropdown(options = family, value = 'Open Sans', disabled = False,
                                   layout = Layout(width='290px', height='25px'))
        
        label_Z = widgets.SelectionSlider(options=range(5, 31),value=12,disabled=False,
                                              description = 'Z label:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        
        ticklabels_Z = widgets.SelectionSlider(options=range(5, 31),value=11,disabled=False,
                                              description = 'Z tick labels:',
                                        continuous_update=False,orientation='horizontal',readout=True,
                                            layout=Layout(width='400px', height='25px'))
        
        
        filename = widgets.Text(value='', placeholder='File name.txt', description='', disabled=False, layout = Layout(width='195px', height='25px'))
        
        filename_plot = widgets.Text(value='', placeholder='Chart name', description='', disabled=False, layout = Layout(width='270px', height='25px'))
        
        cambiar_sample = widgets.ToggleButtons(options=['Code1', 'Code2', 'Code3'], value = 'Code1')
        cambiar_sample.style.button_width = '70px'

        
        def BetA(tipo_kits_1, VariablE, select_PCA):
            
            ElementoS = list(PCA_KITS[tipo_kits_1][VariablE]['category'].keys())
            
            DF_3 = PCA_KITS[tipo_kits_1][VariablE]['df_3']
            sam_clase = dict(zip(DF_3.index, DF_3.clase))
            Permanova = PCA_KITS[tipo_kits_1][VariablE]['permanova']
            
            
            ######

            VAR0 = 0
            VAR1 = 1
            VAR2 = 2
            VAR3 = 3
            VAR4 = 4
            VAR5 = 5
            VAR6 = 6

            width_cuadro = str(100)
            
            blanca3 = widgets.Button(layout=Layout(width='10px', height='2px'), disabled=True)
            blanca3.style.button_color = 'white'
            
            if VariablE == 'Location':
                for unoo in [variables[VAR0]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
            if VariablE == 'ug OTA/kg':
                for unoo in [variables[VAR1]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Cultivation':
                for unoo in [variables[VAR2]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Coffee Variety':
                for unoo in [variables[VAR3]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Genomic DNA kit':
                for unoo in [variables[VAR4]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Drying Time (Days)':
                for unoo in [variables[VAR5]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
                
            if VariablE == 'Postharvest Processing':
                for unoo in [variables[VAR6]]:
                    Lab_Var0 = widgets.Label(unoo)
                    lista_de_colores = list(diccionarios_colores[len(variables_items_ordened[unoo])].keys())
                    SalidA0 = widgets.Dropdown(options = lista_de_colores, value='tab20_v0', layout=Layout(width=width_cuadro+'px', height='28px'))
                    def show_rampa(SalidA0):
                        #barcolor_v2(lista1 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                        patron = SalidA0.split('_v')[0]
                        barcolor_v1(lista1 = colores[patron], lista2 = diccionarios_colores[len(variables_items_ordened[unoo])][SalidA0])
                    OUTshow_rampa0 = widgets.interactive_output(show_rampa, {'SalidA0':SalidA0})
                RAMPA_ELEGIDA = HBox([VBox([blanca3, VBox([Lab_Var0, SalidA0])]), OUTshow_rampa0])
            
            ivl = PCA_KITS[tipo_kits_1][VariablE]['ivl']
            ivw = PCA_KITS[tipo_kits_1][VariablE]['ivw']
            dvw = PCA_KITS[tipo_kits_1][VariablE]['dvw']
            iv_ticks = PCA_KITS[tipo_kits_1][VariablE]['iv_ticks']
            xlines = PCA_KITS[tipo_kits_1][VariablE]['xlines']
            ylines = PCA_KITS[tipo_kits_1][VariablE]['ylines']
            color_list = PCA_KITS[tipo_kits_1][VariablE]['color_list']
            above_threshold_color = PCA_KITS[tipo_kits_1][VariablE]['above_threshold_color']
            id_samples = KITS['Both kits']
            
            
            
            if tipo_kits_1 in ['DNPowerSoil', 'DNMicrobial'] and VariablE == 'Genomic DNA kit':
                print('!!! All values in the grouping vector are the same. !!!')
                
            else:
            
                if select_PCA in ('2D', '3D'):
                    import datetime
                    
                    def BetA2(markers_point, size_point, alfa_point, alfa_linea, ver_num_muestra, ver_dendograma, ver_var_name,
                              family_den, label_X, ticklabels_X, label_Y, ticklabels_Y, family_axis, label_Z, ticklabels_Z, cambiar_sample, SalidA0):
                        
                        category = dict(zip(ElementoS, diccionarios_colores[len(variables_items_ordened[VariablE])][SalidA0]))
                        col = DF_3['clase'].map(category)
                 
                        if select_PCA == '2D':

                            LETRA_2pca = 10
                            plot_alto_pca2 = 6

                            mpl.rcParams.update(mpl.rcParamsDefault)

                            fig = plt.figure(figsize = (plot_alto_pca2, plot_alto_pca2))

                            ax = fig.add_axes([0, 0, 1, 0.8])
                            ax.set_aspect('equal')

                            ax.set_xlim(DF_3.PCA1.min() + (DF_3.PCA1.min() * 0.2), DF_3.PCA1.max() + (DF_3.PCA1.max() * 0.2))
                            ax.set_ylim(DF_3.PCA2.min() + (DF_3.PCA2.min() * 0.2), DF_3.PCA2.max() + (DF_3.PCA2.max() * 0.2))


                            radio = (size_point/5) # el 5 es constante

                            sam_circle_patches= {}
                            sam_circle_patches_pos = {}
                            sam_text_nodos = {}
                            for xx, yy, label, c in zip(DF_3.PCA1, DF_3.PCA2, id_samples, col):

                                if markers_point == '⏺':
                                    CC = ax.add_patch(Circle((xx, yy), radio, angle= 0, facecolor= c, zorder=2, alpha = 0.75, edgecolor = 'white', linewidth = 0.75,
                                                            url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[label]))
                                    sam_circle_patches[label] = CC

                                if markers_point == '⏹':
                                    CC = ax.add_patch(Rectangle((xx- (radio/2), yy - (radio/2)), radio, radio, angle = 0, facecolor=c, zorder=2, alpha = 0.75, edgecolor = 'white', linewidth=0.75,
                                                          url = 'https://www.ncbi.nlm.nih.gov/biosample/?term='+sample_BioSample[label]))
                                    sam_circle_patches[label] = CC


                                AA = ax.text(xx + (radio/2), yy, ' '+label.split('_')[-1], ha = 'left', va = 'center', color = 'black', fontsize = 8, alpha = 1, zorder = 2, family = '3ds Light')
                                sam_text_nodos[label] = AA
                                #ax.scatter(xx, yy, c= c, s = 150, marker = 'o', alpha = 0.8,  edgecolors='white', linewidths = 0.75, zorder=2)

                                sam_circle_patches_pos[label] = (xx, yy)


                            plt.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, zorder = 0)
                            plt.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, zorder = 0)

                            plt.xlabel('PCA 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)', size=LETRA_2pca, weight="bold")
                            plt.ylabel('PCA 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)', size=LETRA_2pca, weight="bold")


                            var_texto = {}
                            lines_plot = {}
                            
                            for cl in DF_3.clase.drop_duplicates():
                                points = []
                                a = DF_3[DF_3.clase == cl].PCA1.values
                                b = DF_3[DF_3.clase == cl].PCA2.values
                                #ax.scatter(np.mean(a), np.mean(b), zorder=0, c= 'white', edgecolors = dict_variable_element_colors[variable][cl], s = 3000, marker = 'o', alpha = 0.1, lw = 5)
                                GG = ax.text(np.mean(a), np.mean(b), cl, ha = 'center', va = 'center', color = 'black', fontsize = 15, alpha = 0.1, weight = 'bold', zorder = 2)
                                var_texto[cl] = GG

                                JJ = []
                                for A, B in zip(a, b):
                                    HH = plt.plot([np.mean(a), A], [np.mean(b), B],
                                             linestyle='--', markeredgewidth=0, zorder=0, markersize=0, color=category[cl], linewidth=1, alpha = 0.3)
                                    JJ.append(HH)

                                lines_plot[cl] = JJ

                            ax.text(ax.get_xlim()[0], ax.get_ylim()[1] + (ax.get_ylim()[1]*0.01), 'PERMANOVA: p='+str(Permanova)+'\n',
                                     fontsize=9, ha='left', va = 'center')


                            proxies = []
                            labels = []
                            for cat in category:
                                proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                                         c=category[cat], marker='o', alpha = 0.75,
                                                         markersize = 7,
                                                         markeredgecolor=category[cat], linewidth = 0)
                                proxies.append(proxy)
                                labels.append(cat)


                            ax.legend(proxies, list(category.keys()), title = ' '.join(["$\\bf{"+i+"}$" for i in VariablE.split(' ')]), numpoints=1, loc=2,borderpad = 0.2,
                                      handletextpad=-0.2,prop={'size':9},
                                              bbox_to_anchor=(1.005, 1.015))


                            ax.tick_params(bottom=True, right=False, top=False, left=True, width = 2, length=4, color='gainsboro')
                            #plt.gca().tick_params(which='major', width = 2, length=4, color='red')
                            for lados in ['left', 'right', 'top', 'bottom']:
                                plt.gca().spines[lados].set_linewidth(2)
                                plt.gca().spines[lados].set_color('gainsboro')


                            ax0 = fig.add_axes([1.3, 0, 0.3, 0.75]) 

                            ax0.set_xlim([dvw, -0.2])
                            ax0.set_ylim([0, ivw])

                            ax0.set_yticks(iv_ticks)

                            ax0.yaxis.set_ticks_position('right')


                            for line in ax0.get_yticklines():
                                line.set_visible(False)

                            colors_used = _remove_dups(color_list)
                            color_to_lines = {}
                            for color in colors_used:
                                color_to_lines[color] = []
                            for (xline, yline, color) in zip(xlines, ylines, color_list):
                                color_to_lines[color].append(list(zip(xline, yline)))

                            colors_to_collections = {}
                            for color in colors_used:
                                coll = matplotlib.collections.LineCollection(color_to_lines[color], colors=(color,))
                                colors_to_collections[color] = coll

                            for color in colors_used:
                                if color != above_threshold_color:
                                    ax0.add_collection(colors_to_collections[color])

                            if above_threshold_color in colors_to_collections:
                                    ax0.add_collection(colors_to_collections[above_threshold_color])

                            for e, i in enumerate(ax0.collections):
                                i.set_color(Set1[e])

                            ax0.text(0, ivw + (ivw*0.05), 'Clustering', fontsize=10, ha='left', va = 'center', weight = 'bold')

                            scatter_ax0 = {}
                            texto_ax0 = {}
                            for i, j in zip(ax0.get_yticks(), ivl):
                                OO = ax0.text(-0.13, i, j, ha = 'left', va = 'center', color = 'black', fontsize = 9, alpha = 1, zorder = 2)
                                texto_ax0[j] = OO
                                LL = ax0.scatter(-0.07, i, c= category[sam_clase[j]], s = 50, marker = 'o', alpha = 1,  edgecolors='none', linewidths = 0.75, zorder=2)
                                scatter_ax0[j] = LL

                            ax0.axis('off')

                            plt.close()


                            def set_locus():

                                if markers_point == '⏺':
                                    [sam_circle_patches[i].set_radius(size_point/5) for i in sam_circle_patches]

                                    proxies = []
                                    labels = []
                                    for cat in category:
                                        proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                                                 c=category[cat], marker='o', alpha = 0.75,
                                                                 markersize = 7,
                                                                 markeredgecolor=category[cat], linewidth = 0)
                                        proxies.append(proxy)
                                        labels.append(cat)


                                    ax.legend(proxies, list(category.keys()), title = ' '.join(["$\\bf{"+i+"}$" for i in VariablE.split(' ')]), numpoints=1, loc=2,borderpad = 0.2,
                                              handletextpad=-0.2,prop={'size':9, 'family':family_axis},
                                                      bbox_to_anchor=(1.005, 1.015))

                                if markers_point == '⏹':
                                    # cambia el tamano de los markers, pero al cambiar el tamano se camvia la posocion
                                    [sam_circle_patches[i].set_width((size_point/5)*1.7) for i in sam_circle_patches]
                                    [sam_circle_patches[i].set_height((size_point/5)*1.7) for i in sam_circle_patches]
                                    # con esto se reajustan las posiciones de los nuevos markers
                                    [sam_circle_patches[i].set_x(sam_circle_patches_pos[i][0] - ((size_point/5)/2)*1.7) for i in sam_circle_patches]
                                    [sam_circle_patches[i].set_y(sam_circle_patches_pos[i][1] - ((size_point/5)/2)*1.7) for i in sam_circle_patches]


                                    proxies = []
                                    labels = []
                                    for cat in category:
                                        proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                                                 c=category[cat], marker='s', alpha = 0.75,
                                                                 markersize = 7,
                                                                 markeredgecolor=category[cat], linewidth = 0)
                                        proxies.append(proxy)
                                        labels.append(cat)


                                    ax.legend(proxies, list(category.keys()), title = ' '.join(["$\\bf{"+i+"}$" for i in VariablE.split(' ')]), title_fontsize = 10, numpoints=1, loc=2,borderpad = 0.2,
                                              handletextpad=-0.2,prop={'size':9, 'family':family_axis},
                                                      bbox_to_anchor=(1.005, 1.015))


                                if ver_dendograma == 'False':
                                    #[scatter_ax0[i].set_alpha(0) for i in scatter_ax0]
                                    #[texto_ax0[i].set_alpha(0) for i in texto_ax0]
                                    #for e, i in enumerate(ax0.collections):
                                    #    i.set_linewidth(0)
                                    ax0.clear()
                                    ax0.axis('off')

                                if ver_dendograma == 'True':
                                    pass

                                if ver_num_muestra == 'False':
                                    [sam_text_nodos[i].set_alpha(0) for i in sam_text_nodos]
                                if ver_num_muestra == 'True':
                                    pass


                                [scatter_ax0[i].set_paths([marcas[{'⏺':'o', '⏹':'s'}[markers_point]]]) for i in scatter_ax0]

                                [var_texto[i].set_alpha(ver_var_name) for i in var_texto]


                                [sam_circle_patches[i].set_alpha(alfa_point) for i in sam_circle_patches]

                                [j[0].set_alpha(alfa_linea) for i in lines_plot for j in lines_plot[i]]
                                
                                [texto_ax0[i].set_family(family_den) for i in texto_ax0]
                                
                                for labejex in ax.xaxis.get_ticklabels():
                                    labejex.set_fontsize(ticklabels_X)
                                for labejey in ax.yaxis.get_ticklabels():
                                    labejey.set_fontsize(ticklabels_Y)
                                for tickx in ax.xaxis.get_major_ticks():
                                    tickx.label.set_fontfamily([family_axis]) 
                                for ticky in ax.yaxis.get_major_ticks():
                                    ticky.label.set_fontfamily([family_axis])
                                ax.xaxis.get_label().set_fontsize(label_X)
                                ax.yaxis.get_label().set_fontsize(label_Y)
                                ax.xaxis.get_label().set_fontfamily([family_axis])
                                ax.yaxis.get_label().set_fontfamily([family_axis])
                                
                                if cambiar_sample == 'Code1':
                                    pass
                                if cambiar_sample == 'Code2':
                                    [texto_ax0[i].set_text(name_code[texto_ax0[i].get_text()]) for i in texto_ax0]
                                if cambiar_sample == 'Code3':
                                    [texto_ax0[i].set_text(name_code2[texto_ax0[i].get_text()]) for i in texto_ax0]


                                display(fig)

                        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                        if select_PCA == '3D':

                            plot_alto_pca3 = 9
                            LETRA_3pca = 9

                            mpl.rcParams.update(mpl.rcParamsDefault)
                            plt.rcParams['grid.color'] = "whitesmoke"


                            fig = plt.figure(figsize=(plot_alto_pca3, plot_alto_pca3-2))

                            ax = fig.add_subplot(111, projection = '3d', facecolor = 'white')

                            ax.grid(True)
                            #ax.xaxis.pane.fill = False
                            #ax.yaxis.pane.fill = False
                            #ax.zaxis.pane.fill = False
                            
                            sam_text_nodos = {}
                            sam_circle_patches = {}
                            for x, y, z, label, c  in zip(DF_3.PCA1, DF_3.PCA2, DF_3.PCA3, id_samples, col):
                                FF = ax.scatter(x, y, z, c= c, s = 100, marker = 'o', alpha = 0.7, edgecolors='white', linewidths = 1, zorder=1)
                                sam_circle_patches[label] = FF
                                
                                
                                AA = ax.text(x, y, z, '   '+label.split('_')[-1], ha = 'left', va = 'center', color = 'black', fontsize = 9, alpha = 1, zorder = 2, family = '3ds Light')
                                sam_text_nodos[label] = AA
                                

                            ax.text(ax.get_xlim()[1] * 0.1, ax.get_ylim()[1], ax.get_zlim()[1], 'PERMANOVA: p='+str(Permanova)+'\n\n\n\n',
                                    fontsize=9, ha='left', va = 'center')

                            proxies = []
                            labels = []
                            for cat in category:
                                proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                                         c=category[cat], marker={'⏺':'o', '⏹':'s'}[markers_point], alpha = 0.7,
                                                         markersize = 7,
                                                         markeredgecolor=category[cat], linewidth = 0)
                                proxies.append(proxy)
                                labels.append(cat)


                            ax.legend(proxies, list(category.keys()), title = ' '.join(["$\\bf{"+i+"}$" for i in variable.split(' ')]), numpoints=1, loc=2,borderpad = 0.2,
                                      handletextpad=-0.2,prop={'size':9, 'family':family_axis},
                                              bbox_to_anchor=(1.015, 0.93))
                            
                            var_texto = {}
                            lines_plot = {}
                            for cl in DF_3.clase.drop_duplicates():
                                points = []
                                a = DF_3[DF_3.clase == cl].PCA1.values
                                b = DF_3[DF_3.clase == cl].PCA2.values
                                c = DF_3[DF_3.clase == cl].PCA3.values
                                #ax.scatter(np.mean(a), np.mean(b), np.mean(c), zorder=0, c= 'black', s = 5, marker = 'o', alpha = 1, lw = 0)
                                GG = ax.text(np.mean(a), np.mean(b), np.mean(c), cl, ha = 'center', va = 'center', color = 'black', fontsize = 15, alpha = 0.1, weight = 'bold', zorder = 2)
                                var_texto[cl] = GG
                                
                                JJ = []
                                for A, B, C in zip(a, b, c):
                                    HH = plt.plot([np.mean(a), A], [np.mean(b), B], [np.mean(c), C],
                                             linestyle='--', markeredgewidth=0, zorder=0, markersize=0, color=category[cl], linewidth=0.5, alpha = 0.3)
                                    JJ.append(HH)
                                lines_plot[cl] = JJ

                            ax.set_xlabel('PCA 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)', size=LETRA_3pca, weight="bold")
                            ax.set_ylabel('PCA 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)', size=LETRA_3pca, weight="bold")
                            ax.set_zlabel('PCA 3 ('+str(round(pca3.explained_variance_ratio_[2]*100,2))+'%)', size=LETRA_3pca, weight="bold")

                            plt.xticks(size=9)
                            plt.yticks(size=9)

                            ax.zaxis.set_tick_params(labelsize=9)

                            fccolor3D = 'white'
                            ax.xaxis.pane.set_facecolor(fccolor3D)
                            ax.yaxis.pane.set_facecolor(fccolor3D)
                            ax.zaxis.pane.set_facecolor(fccolor3D)
                            ax.xaxis.pane.set_edgecolor('gainsboro')
                            ax.yaxis.pane.set_edgecolor('gainsboro')
                            ax.zaxis.pane.set_edgecolor('gainsboro')
                            ax.xaxis.pane.set_alpha(1)
                            ax.yaxis.pane.set_alpha(1)
                            ax.zaxis.pane.set_alpha(1)
                            ax.xaxis.pane.set_linewidth(2)
                            ax.yaxis.pane.set_linewidth(2)
                            ax.zaxis.pane.set_linewidth(2)

                            #---------------------------------------------

                            ax0 = fig.add_axes([1.18, 0.1, 0.2, 0.65]) 

                            ax0.set_xlim([dvw, -0.2])
                            ax0.set_ylim([0, ivw])

                            ax0.set_yticks(iv_ticks)

                            ax0.yaxis.set_ticks_position('right')


                            for line in ax0.get_yticklines():
                                line.set_visible(False)

                            colors_used = _remove_dups(color_list)
                            color_to_lines = {}
                            for color in colors_used:
                                color_to_lines[color] = []
                            for (xline, yline, color) in zip(xlines, ylines, color_list):
                                color_to_lines[color].append(list(zip(xline, yline)))

                            colors_to_collections = {}
                            for color in colors_used:
                                coll = matplotlib.collections.LineCollection(color_to_lines[color], colors=(color,))
                                colors_to_collections[color] = coll

                            for color in colors_used:
                                if color != above_threshold_color:
                                    ax0.add_collection(colors_to_collections[color])

                            if above_threshold_color in colors_to_collections:
                                    ax0.add_collection(colors_to_collections[above_threshold_color])

                            for e, i in enumerate(ax0.collections):
                                i.set_color(Set1[e])

                            ax0.text(0, ivw + (ivw*0.05), 'Clustering', fontsize=10, ha='left', va = 'center', weight = 'bold')
                            
                            scatter_ax0 = {}
                            texto_ax0 = {}
                            for i, j in zip(ax0.get_yticks(), ivl):
                                OO = ax0.text(-0.13, i, j, ha = 'left', va = 'center', color = 'black', fontsize = 9, alpha = 1, zorder = 2)
                                texto_ax0[j] = OO
                                LL = ax0.scatter(-0.07, i, c= category[sam_clase[j]], s = 50, marker = 'o', alpha = 1,  edgecolors='none', linewidths = 0.75, zorder=2)
                                scatter_ax0[j] = LL
                            
                            ax0.axis('off')

                            plt.close()

                            def set_locus():
                                

                                [j[0].set_alpha(alfa_linea) for i in lines_plot for j in lines_plot[i]]
                                
                                [sam_circle_patches[i].set_paths([marcas[{'⏺':'o', '⏹':'s'}[markers_point]]]) for i in sam_circle_patches]
                                
                                [scatter_ax0[i].set_paths([marcas[{'⏺':'o', '⏹':'s'}[markers_point]]]) for i in scatter_ax0]
                                
                                [texto_ax0[i].set_family(family_den) for i in texto_ax0]
                                
                                [sam_circle_patches[i].set_sizes([(size_point*2)*100]) for i in sam_circle_patches]
                                
                                [sam_circle_patches[i].set_alpha(alfa_point) for i in sam_circle_patches]
                                
                                [var_texto[i].set_alpha(ver_var_name) for i in var_texto]
                                
                                if cambiar_sample == 'Code1':
                                    pass
                                if cambiar_sample == 'Code2':
                                    [texto_ax0[i].set_text(name_code[texto_ax0[i].get_text()]) for i in texto_ax0]
                                if cambiar_sample == 'Code3':
                                    [texto_ax0[i].set_text(name_code2[texto_ax0[i].get_text()]) for i in texto_ax0]
                                
                                if ver_dendograma == 'False':
                                    ax0.clear()
                                    ax0.axis('off')

                                if ver_dendograma == 'True':
                                    pass
                                
                                for labejex in ax.xaxis.get_ticklabels():
                                    labejex.set_fontsize(ticklabels_X)
                                for labejey in ax.yaxis.get_ticklabels():
                                    labejey.set_fontsize(ticklabels_Y)
                                for labejey in ax.zaxis.get_ticklabels():
                                    labejey.set_fontsize(ticklabels_Z)
                                    
                                    
                                for tickx in ax.xaxis.get_major_ticks():
                                    tickx.label.set_fontfamily([family_axis]) 
                                for ticky in ax.yaxis.get_major_ticks():
                                    ticky.label.set_fontfamily([family_axis])
                                for ticky in ax.zaxis.get_major_ticks():
                                    ticky.label.set_fontfamily([family_axis])    
                                    
                                ax.xaxis.get_label().set_fontsize(label_X)
                                ax.yaxis.get_label().set_fontsize(label_Y)
                                ax.zaxis.get_label().set_fontsize(label_Z)
                                
                                ax.xaxis.get_label().set_fontfamily([family_axis])
                                ax.yaxis.get_label().set_fontfamily([family_axis])
                                ax.zaxis.get_label().set_fontfamily([family_axis])
                                
                                if ver_num_muestra == 'False':
                                    [sam_text_nodos[i].set_alpha(0) for i in sam_text_nodos]
                                if ver_num_muestra == 'True':
                                    pass
                                

                                display(fig)
                        
                        
                        
                        # --- save
                        ancho_plot_save = str(71)
                        png1 = widgets.Button(description="PNG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
                        png1.style.button_color = 'gold'
                        output1 = widgets.Output()
                        def button_clicked1(b):
                            with output1:
                                clear_output(True)
                                if filename_plot.value is '':
                                    nombre_grafico = 'Beta_diversity_phyl_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                                if filename_plot.value is not '':
                                    nombre_grafico = filename_plot.value
                                    
                                fig.savefig('plots_asv/beta_diversity/'+nombre_grafico+'.png', dpi = 900, bbox_inches= 'tight')
                                
                                
                        png1.on_click(button_clicked1)
                        #----
                        jpeg1 = widgets.Button(description="JPEG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
                        jpeg1.style.button_color = 'gold'
                        output2 = widgets.Output()
                        def button_clicked2(b):
                            with output2:
                                clear_output(True)
                                if filename_plot.value is '':
                                    nombre_grafico = 'Beta_diversity_phyl_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                                if filename_plot.value is not '':
                                    nombre_grafico = filename_plot.value
                                    
                                fig.savefig('plots_asv/beta_diversity/'+nombre_grafico+'.jpeg', dpi = 900, bbox_inches= 'tight')
                                
                        jpeg1.on_click(button_clicked2)
                        #----
                        svg1 = widgets.Button(description="SVG", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
                        svg1.style.button_color = 'gold'
                        output3 = widgets.Output()
                        def button_clicked3(b):
                            with output3:
                                clear_output(True)
                                if filename_plot.value is '':
                                    nombre_grafico = 'Beta_diversity_phyl_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                                if filename_plot.value is not '':
                                    nombre_grafico = filename_plot.value
                                    
                                fig.savefig('plots_asv/beta_diversity/'+nombre_grafico+'.svg', dpi = 900, bbox_inches= 'tight')
                                
                        svg1.on_click(button_clicked3)
                        #----
                        pdf1 = widgets.Button(description="PDF", icon = 'fa-bar-chart', layout=Layout(width=ancho_plot_save+'px'))
                        pdf1.style.button_color = 'gold'
                        output4 = widgets.Output()
                        def button_clicked4(b):
                            with output4:
                                clear_output(True)
                                if filename_plot.value is '':
                                    nombre_grafico = 'Beta_diversity_phyl_chart_'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')
                                if filename_plot.value is not '':
                                    nombre_grafico = filename_plot.value
                                    
                                fig.savefig('plots_asv/beta_diversity/'+nombre_grafico+'.pdf', dpi = 900, bbox_inches= 'tight')
                                
                        pdf1.on_click(button_clicked4)
                        ### ---


                        params = widgets.Button(description="SAVE", icon = 'fa-save', layout=Layout(width=ancho_plot_save+'px'))
                        params.style.button_color = 'cyan'
                        output5 = widgets.Output()
                        def button_clicked5(b):
                            with output5:
                                clear_output(True)

                                if filename.value is '':
                                    nombre_archivo = 'Beta_diversity_phyl_chart_params_ASVs.txt'
                                if filename.value is not '':
                                    nombre_archivo = filename.value

                                with open('plots_asv/beta_diversity/'+nombre_archivo, 'w') as fq:
                                    fq.write('#Saved parameters\n')
                                    fq.write('Clustering:ASVs\n' )
                                    """
                                    diccionario para guardar en un archivo los parametros usados para volver a reproducirlos posteriormente
                                    """
                                    parametros = {'Kit:':tipo_kits_1,
                                                  'Variable:':VariablE,
                                                  'Chart type:':select_PCA,
                                                  'Marker:':{'⏺':'circle', '⏹':'square'}[markers_point],
                                                  'Sample number:':ver_num_muestra,
                                                  'Marker size:':size_point,
                                                  'Marked alpha:':alfa_point,
                                                  'Line alpha:':alfa_linea,
                                                  'Text alpha:':ver_var_name,
                                                  'Dendrogram font:':family_den,
                                                  'Dendrogram:':ver_dendograma,
                                                  'Rename sample:':cambiar_sample,
                                                  'X label:':label_X,
                                                  'X tick label:':ticklabels_X,
                                                  'Y label:':label_Y,
                                                  'Y tick labels:':ticklabels_Y,
                                                  'Z label:':label_Y,
                                                  'Z tick labels:':ticklabels_Y,
                                                  'Axis font:':family_axis,
                                                  '\n###\nWarning:':'The variables take the same colors assigned in the taxonomy analysis.'}
                                    for w in parametros:
                                        fq.write(w+str(parametros[w])+'\n')

                        params.on_click(button_clicked5)
                        #------
                        
                        items_save_1 = VBox([HBox([widgets.HTML('<font color = grey> <b style="font-size:0.7vw">SAVE CHART: </b>'), blanca, filename_plot]),
                                             HBox([blanca, widgets.Label('Formats:'), png1, jpeg1, svg1, pdf1]),
                                     HBox([blanca, widgets.Label('Chart parameters:'), filename, params])])
                        items_save_1_box = Box(children=[items_save_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_save_1.children) * 34))+'px'))
                        
                        display(VBox([HBox([RAMPA_ELEGIDA, items_save_1_box]), output5]))
                        
                        set_locus()
                        

                    OUT_BetA2 = widgets.interactive_output(BetA2, {'markers_point':markers_point, 'size_point':size_point, 'ver_num_muestra':ver_num_muestra,
                                                     'ver_dendograma':ver_dendograma, 'ver_var_name':ver_var_name,
                                                     'alfa_point':alfa_point, 'alfa_linea':alfa_linea,
                                                                  'label_X':label_X, 'ticklabels_X':ticklabels_X, 'label_Y':label_Y, 'ticklabels_Y':ticklabels_Y, 'family_axis':family_axis,
                                                                   'label_Z':label_Z, 'ticklabels_Z':ticklabels_Z, 'family_den':family_den, 'cambiar_sample':cambiar_sample, 'SalidA0':SalidA0})
                    
                    
                    display(OUT_BetA2)


                if select_PCA == '3Di':
                    import datetime
                    category = dict_variable_element_colors[VariablE]
                    
                    blanca0 = widgets.Button(layout=Layout(width='3px', height='25px'), disabled=True)
                    blanca0.style.button_color = 'white'

                    # --- save
                    html1 = widgets.Button(description="HTML", icon = 'fa-bar-chart', layout=Layout(width='71px'))
                    html1.style.button_color = 'gold'
                    output11 = widgets.Output()
                    def button_clicked11(b):
                        with output11:
                            clear_output(True)
                            plotly.offline.plot(FigurA, filename = 'plots_asv/beta_diversity/'+tipo_kits_1+'_'+re.sub('[/]|[(]|[)]|[*]|[|]', '_', VariablE)+'_ASVs'+datetime.datetime.now().strftime('%d.%B.%Y_%I-%M%p')+'.html', auto_open=False)

                    html1.on_click(button_clicked11)

                    ver_linea = widgets.ToggleButtons(options=['True', 'False'], value = 'True')
                    ver_linea.style.button_width = '55px'
                    ver_texto = widgets.ToggleButtons(options=['True', 'False'], value = 'True')
                    ver_texto.style.button_width = '55px'
                    ver_estilo = widgets.ToggleButtons(options=['White', 'Black'], value = 'White')
                    ver_estilo.style.button_width = '55px'


                    desc_ver_linea = HBox([widgets.Label('Show line:'), ver_linea, blanca0, widgets.Label('Sample number:'), ver_texto, blanca0, widgets.Label('Style:'), ver_estilo, blanca0, widgets.Label('Save 3Di:'), html1])
                    plot_3di = HBox([widgets.HTML('<font color = black> <b style="font-size:0.6vw">3Di Settings: </b>'), blanca0, desc_ver_linea])
                    set_box = Box(children=[plot_3di], layout=Layout(border='1px solid black', width='898px', height='37px'))



                    Lines_poS = {}
                    for cl in DF_3.clase.drop_duplicates():
                        a = DF_3[DF_3.clase == cl].PCA1.values
                        b = DF_3[DF_3.clase == cl].PCA2.values
                        c = DF_3[DF_3.clase == cl].PCA3.values

                        lines_positions_data = {}
                        n = 0
                        for A, B, C in zip(a, b, c):
                            n += 1        
                            lines_positions_data['data'+str(n)] = {'x':np.array([np.mean(a), A]), 'y':np.array([np.mean(b), B]), 'z':np.array([np.mean(c), C])}
                        Lines_poS[cl] = lines_positions_data


                    DF3D = DF_3
                    DF3D['Size'] = 15
                    DF3D[VariablE] = DF3D.clase
                    DF3D['Sample'] = DF3D.index
                    DF3D = DF3D.reset_index(drop = True)

                    colorr = list(category.values())

                    labels = {'PCA1': 'PC 1 ('+str(round(pca3.explained_variance_ratio_[0]*100,2))+'%)',
                              'PCA2': 'PC 2 ('+str(round(pca3.explained_variance_ratio_[1]*100,2))+'%)',
                              'PCA3': 'PC 3 ('+str(round(pca3.explained_variance_ratio_[2]*100,2))+'%)'}

                    fig = go.Figure()

                    scatter_clases = 0
                    for i in DF3D.clase.unique():
                        x = DF3D[DF3D.clase == i].PCA1.values
                        y = DF3D[DF3D.clase == i].PCA2.values
                        z = DF3D[DF3D.clase == i].PCA3.values
                        lab = DF3D[DF3D.clase == i].Sample.tolist()
                        fig.add_trace(go.Scatter3d(x=x, y=y, z=z,
                                                           name = i,
                                                           projection=None,
                                                           text= 'none',
                                                           hovertext = [VariablE+'='+i+'<br>'+labels['PCA1']+'='+str(round(a, 5))+'<br>'+labels['PCA2']+'='+str(round(b, 5))+'<br>'+labels['PCA3']+'='+str(round(c, 5))+'<br>'+k for k, a, b, c in zip(lab, x, y, z)],
                                                           hoverinfo='text',
                                                           marker_color = category[i],
                                                           opacity=1,
                                                           showlegend=True,
                                                           marker=dict(size=6, opacity=1, line=dict(color=category[i], width=200)),
                                                           textfont=dict(family="3ds Light", size=10, color=category[i]),
                                                           mode='markers'))
                        scatter_clases += 1


                    scatter_lines = scatter_clases
                    for h in Lines_poS:
                        for i in Lines_poS[h]:
                            fig.add_trace(go.Scatter3d(x=Lines_poS[h][i]['x'], y=Lines_poS[h][i]['y'], z=Lines_poS[h][i]['z'], hoverinfo='none', showlegend=False,
                                                        opacity=0.5, mode='lines', name='lines', line = dict(color=category[h], width=5, dash='solid')))
                            scatter_lines += 1


                    scater_con_text = scatter_lines
                    for i in DF3D.clase.unique():
                        x = DF3D[DF3D.clase == i].PCA1.values
                        y = DF3D[DF3D.clase == i].PCA2.values
                        z = DF3D[DF3D.clase == i].PCA3.values
                        lab = DF3D[DF3D.clase == i].Sample.tolist()
                        fig.add_trace(go.Scatter3d(x=x, y=y, z=z,
                                                           name = i,
                                                           text = [i.split('_')[-1] for i in lab],
                                                           projection=None,
                                                           hovertext = [VariablE+'='+i+'<br>'+labels['PCA1']+'='+str(round(a, 5))+'<br>'+labels['PCA2']+'='+str(round(b, 5))+'<br>'+labels['PCA3']+'='+str(round(c, 5))+'<br>'+k for k, a, b, c in zip(lab, x, y, z)],
                                                           hoverinfo='text',
                                                           marker_color = category[i],
                                                           opacity=1,
                                                           showlegend=False,
                                                           marker=dict(size=6, opacity=1, line=dict(color=category[i], width=200)),
                                                           textfont=dict(family="3ds Light", size=13),
                                                           mode='markers+text'))
                        scater_con_text += 1


                    fig.update_layout(legend_title_text=VariablE, autosize=True,width=900, height=650,margin=dict(l=0, r=5, b=0, t=0, pad=0),
                                      scene = dict(xaxis_title=labels['PCA1'], yaxis_title=labels['PCA2'], zaxis_title=labels['PCA3']),
                                         font_size=13, template="plotly_white", legend=dict(yanchor="top", y=0.9, xanchor="left", x=1))

                    FigurA = go.FigureWidget(fig)
                    display(VBox([set_box, FigurA]))

                    def res(f):
                        with FigurA.batch_update():
                            if ver_linea.value == 'True':
                                for i in FigurA.data[scatter_clases:scatter_lines]:
                                    i.opacity = 0.5
                            if ver_linea.value == 'False':
                                for i in FigurA.data[scatter_clases:scatter_lines]:
                                    i.opacity = 0

                            if ver_texto.value == 'True':
                                for i in FigurA.data[scatter_lines:]:
                                    i.mode = 'markers+text'
                            if ver_texto.value == 'False':
                                for i in FigurA.data[scatter_lines:]:
                                    i.mode = 'markers'
                                    
                            if ver_estilo.value == 'White':
                                FigurA.layout.template = "plotly_white"
                                
                            if ver_estilo.value == 'Black':
                                FigurA.layout.template = "plotly_dark"


                    ver_linea.observe(res, names="value")
                    ver_texto.observe(res, names="value")
                    ver_estilo.observe(res, names="value")
                
        OUT_BetA = widgets.interactive_output(BetA, {'tipo_kits_1':tipo_kits_1, 'VariablE':VariablE, 'select_PCA':select_PCA,})
        
        items_data = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'DATA']
        items_axis = [widgets.HTML('<font color = grey> <b style="font-size:0.7vw">'+i+'</b>') for i in 'AXIS']
        items_save = [widgets.HTML('<font color = white> <b style="font-size:0.7vw">'+i+'</b>') for i in 'S']
        
        
        
        
        items_plot = VBox([HBox([blanca, widgets.Label('Marker size:'),
                                 markers_point,
                                 widgets.Label('Sample number:'), ver_num_muestra]),
                                 size_point,
                                 alfa_point,
                                 alfa_linea,
                                 ver_var_name,
                                 HBox([widgets.Label('Dendrogram font:'), family_den]),
                                 HBox([widgets.Label('Dendrogram:'), ver_dendograma]),
                                 HBox([widgets.Label('Rename sample:'), cambiar_sample])])
        
        items_plot_box = Box(children=[items_plot], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_plot.children) * 31))+'px'))
        
        items_axis_1 = VBox([label_X,
                                 ticklabels_X,
                                 label_Y,
                                 ticklabels_Y,
                                 label_Z,
                                 ticklabels_Z,
                                 HBox([blanca, widgets.Label('Axis font:'), family_axis])])
        
        items_axis_1_box = Box(children=[items_axis_1], layout=Layout(border='1px solid gainsboro', width='410px', height=str(int(len(items_axis_1.children) * 31))+'px'))
        
        
        display(HBox([VBox([VBox([widgets.HTML('<font color = #1976d2> <b style="font-size:0.8vw">KITS</b>'), tipo_kit_1_box]),
                            VBox([widgets.HTML('<font color = gray> <b style="font-size:0.8vw">VARIABLES</b>'), tipo_variable_box]),
                            VBox([widgets.HTML('<font color = gray> <b style="font-size:0.6vw">Chart type:</b>'), boton_select_PCA])]), blanca,
                      VBox([widgets.HTML('<font color = grey> <i class="fa fa-cog fa-2x fa-fw"></i> <b style="font-size:0.8vw">PLOT SETTINGS (2D & 3D)</b>'),
                            HBox([Box(children =[VBox(items_data)]), items_plot_box]),
                             
                            HBox([Box(children =[VBox(items_axis)]), items_axis_1_box]),
                            
                           
                           ]),
                      blanca, OUT_BetA]))
        
BETA_PHY_button.on_click(button_clicked)





BETA_DIVERSITY_PHYL = VBox([HBox([BETA_PHY_button, estatico_fil_box]), BETA_PHY_output])























graficas = {'Rarefaction Analysis': RARE_ANALYSIS,
            'Alpha Diversity': DIVERSITY_CALCULATIONS,
            'Richness': RICHNESS_CALCULATIONS,
            'Taxonomy': TAXONOMY_ANALYSIS,
            'Beta Diversity (No phylogenetic)': BETA_DIVERSITY_NO_PHYL,
            'Beta Diversity (Phylogenetic)': BETA_DIVERSITY_PHYL}
exe = widgets.ToggleButtons(options=list(graficas.keys()),disabled=False,button_style='warning')
exe.style.button_width = '258px'
exe.style.font_weight = 'bold'





box_exe = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='5px solid gainsboro',
                    width='1700px',
                   height='45px')



def ss(exe):
    display(graficas[exe])
OU = widgets.interactive_output(ss, {'exe':exe})


RESULTS_16S = VBox([Box(children = [HBox([widgets.HTML('<b style="font-size:0.9vw">ANALYSIS: </b>'), exe])], layout = box_exe), OU])





threshold_box = Box(children=[VBox([blanca2, threshold])], layout=Layout(border='5px solid beige', width='1700px', height='150px'))





METAGENOMIC_16S_ASVs = VBox([threshold_box, blanca2, RESULTS_16S])
