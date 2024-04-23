#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
20.Created on Mon Feb 28 15:00:01 2022

@author: student
"""
import colorsys
from turtle import color
from Bio.PDB.PDBList import PDBList
from Bio.SeqUtils import seq1
import mdtraj as mdt
import os
import numpy as np
from random import randrange
import os
import glob
import shutil
import re
tool_files = os.getcwd()

import time
import threading
import csv
from PyQt5 import QtCore, QtWidgets, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QTableWidgetItem
from PyQt5.QtWidgets import QFileDialog ,QMessageBox
from PyQt5.QtWidgets import *
import sys
from colorama import  Fore, Back, Style
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure



import random

global t
global dir_name1
global b
b = ''
global c
c = ''
global x
x = ''
global idscore
idscore = ''
global found
from seq_retrieve import StructSeqRetrieve
try:
    import httplib
except:
    import http.client as httplib





from read_write_seqs import extract_seq
from read_write_seqs import write_seq

import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np

import itertools
import pandas as pd

from apply_matrix import open_file1
#from plot_file import plots
from Bio import AlignIO

import xlwt
from xlwt import Workbook

import webbrowser

import urllib.request
class Ui_V_Motif(object):
    def setupUi(self, V_Motif):
        V_Motif.setObjectName("V_Motif")
        V_Motif.resize(941, 670)
        V_Motif.setMinimumSize(QtCore.QSize(941, 670))
        V_Motif.setMaximumSize(QtCore.QSize(941, 670))
        font = QtGui.QFont()
        font.setFamily("Z003 [urw]")
        font.setPointSize(25)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        V_Motif.setFont(font)
        V_Motif.setMouseTracking(False)
        V_Motif.setTabletTracking(False)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("Icon.jpg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        V_Motif.setWindowIcon(icon)
        V_Motif.setStyleSheet("background-color: rgb(102, 138, 139);")
        V_Motif.setIconSize(QtCore.QSize(50, 50))
        self.centralwidget = QtWidgets.QWidget(V_Motif)
        self.centralwidget.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.centralwidget.setObjectName("centralwidget")
        self.stackedWidget = QtWidgets.QStackedWidget(self.centralwidget)
        self.stackedWidget.setGeometry(QtCore.QRect(-40, -20, 1031, 701))
        self.stackedWidget.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.stackedWidget.setObjectName("stackedWidget")
        self.Page_1 = QtWidgets.QWidget()
        self.Page_1.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.Page_1.setObjectName("Page_1")
        self.Tabs = QtWidgets.QTabWidget(self.Page_1)
        self.Tabs.setGeometry(QtCore.QRect(40, 20, 941, 701))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(15)
        font.setBold(False)
        font.setItalic(True)
        font.setWeight(50)
        self.Tabs.setFont(font)
        self.Tabs.setFocusPolicy(QtCore.Qt.TabFocus)
        self.Tabs.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.Tabs.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.Tabs.setElideMode(QtCore.Qt.ElideMiddle)
        self.Tabs.setTabsClosable(False)
        self.Tabs.setTabBarAutoHide(False)
        self.Tabs.setObjectName("Tabs")
        self.Tab1 = QtWidgets.QWidget()
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.Tab1.setFont(font)
        self.Tab1.setContextMenuPolicy(QtCore.Qt.NoContextMenu)
        self.Tab1.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.Tab1.setObjectName("Tab1")
        self.VMotif = QtWidgets.QLabel(self.Tab1)
        self.VMotif.setGeometry(QtCore.QRect(240, 180, 361, 91))
        font = QtGui.QFont()
        font.setFamily("Z003 [urw]")
        font.setPointSize(50)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(65)
        self.VMotif.setFont(font)
        self.VMotif.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.VMotif.setAlignment(QtCore.Qt.AlignCenter)
        self.VMotif.setObjectName("VMotif")
        self.VMotif_2 = QtWidgets.QLabel(self.Tab1)
        self.VMotif_2.setGeometry(QtCore.QRect(10, 290, 931, 71))
        font = QtGui.QFont()
        font.setFamily("Z003 [urw]")
        font.setPointSize(40)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(45)
        self.VMotif_2.setFont(font)
        self.VMotif_2.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.VMotif_2.setAlignment(QtCore.Qt.AlignCenter)
        self.VMotif_2.setObjectName("VMotif_2")
        self.Tabs.addTab(self.Tab1, "")
        self.Tab2 = QtWidgets.QWidget()
        self.Tab2.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.Tab2.setObjectName("Tab2")
        self.L_Protein_ID = QtWidgets.QLabel(self.Tab2)
        self.L_Protein_ID.setGeometry(QtCore.QRect(230, 150, 201, 41))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(17)
        font.setBold(False)
        font.setItalic(True)
        font.setWeight(50)
        self.L_Protein_ID.setFont(font)
        self.L_Protein_ID.setObjectName("L_Protein_ID")
        self.Enter_ID = QtWidgets.QTextEdit(self.Tab2)
        self.Enter_ID.setGeometry(QtCore.QRect(450, 150, 151, 51))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(15)
        font.setItalic(True)
        self.Enter_ID.setFont(font)
        self.Enter_ID.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.Enter_ID.setObjectName("Enter_ID")
        self.Enter_Output = QtWidgets.QLineEdit(self.Tab2)
        self.Enter_Output.setGeometry(QtCore.QRect(450, 230, 291, 51))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(15)
        font.setItalic(True)
        self.Enter_Output.setFont(font)
        self.Enter_Output.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.Enter_Output.setObjectName("Enter_Output")
        self.B_Execute = QtWidgets.QPushButton(self.Tab2)
        self.B_Execute.setGeometry(QtCore.QRect(410, 320, 250, 51))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(20)
        font.setItalic(True)
        self.B_Execute.setFont(font)
        self.B_Execute.setMouseTracking(False)
        self.B_Execute.setTabletTracking(False)
        self.B_Execute.setFocusPolicy(QtCore.Qt.NoFocus)
        self.B_Execute.setStyleSheet("background-color: rgb(204, 204, 204);\n"
"border-color: rgb(27, 27, 27);")
        self.B_Execute.setObjectName("B_Execute")
        
        self.B_Execute.clicked.connect(self.page_change)
        
        self.L_Output_F = QtWidgets.QPushButton(self.Tab2)
        self.L_Output_F.setGeometry(QtCore.QRect(190, 235, 250, 41))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.L_Output_F.sizePolicy().hasHeightForWidth())
        self.L_Output_F.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(17)
        font.setItalic(True)
        self.L_Output_F.setFont(font)
        self.L_Output_F.setCheckable(False)
        self.L_Output_F.setObjectName("L_Output_F")
        
        self.L_Output_F.clicked.connect(self.get_out)
        
        self.Tabs.addTab(self.Tab2, "")
        
        self.Tab3 = QtWidgets.QWidget()
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.Tab3.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.Tab3.setObjectName("Tabe3")
        self.ProDom = QtWidgets.QCommandLinkButton(self.Tab3)
        self.ProDom.setGeometry(QtCore.QRect(160, 345, 221, 61))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(20)
        font.setBold(False)
        font.setItalic(True)
        font.setWeight(50)
        self.ProDom.setFont(font)
        self.ProDom.setObjectName("ProDom")
        
        self.ProDom.clicked.connect(self.Pro_Dom)
        
        self.MEME = QtWidgets.QCommandLinkButton(self.Tab3)
        self.MEME.setGeometry(QtCore.QRect(660, 345, 211, 61))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(20)
        font.setBold(False)
        font.setItalic(True)
        font.setWeight(50)
        self.MEME.setFont(font)
        self.MEME.setObjectName("MEME")
        
        self.MEME.clicked.connect(self.ME_ME)
        
        self.PRINTS = QtWidgets.QCommandLinkButton(self.Tab3)
        self.PRINTS.setGeometry(QtCore.QRect(420, 345, 221, 61))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(20)
        font.setBold(False)
        font.setItalic(True)
        font.setWeight(50)
        self.PRINTS.setFont(font)
        self.PRINTS.setObjectName("PRINTS")
        
        self.PRINTS.clicked.connect(self.PRI_NTS)
        
        self.Databases = QtWidgets.QLabel(self.Tab3)
        self.Databases.setGeometry(QtCore.QRect(210, 155, 621, 51))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(22)
        font.setItalic(True)
        self.Databases.setFont(font)
        self.Databases.setAlignment(QtCore.Qt.AlignCenter)
        self.Databases.setObjectName("Databases")
        self.Databases_2 = QtWidgets.QLabel(self.Tab3)
        self.Databases_2.setGeometry(QtCore.QRect(210, 230, 621, 51))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(22)
        font.setItalic(True)
        self.Databases_2.setFont(font)
        self.Databases_2.setAlignment(QtCore.Qt.AlignCenter)
        self.Databases_2.setObjectName("Databases_2")

        self.Tabs.addTab(self.Tab3, "")
        
        self.Tab4 = QtWidgets.QWidget()
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.Tab4.setFont(font)
        self.Tab4.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.Tab4.setObjectName("Tab4")
        self.Pic = QtWidgets.QLabel(self.Tab4)
        self.Pic.setGeometry(QtCore.QRect(170, 20, 600, 456))
        self.Pic.setText("")
        self.Pic.setPixmap(QtGui.QPixmap("Icon.jpg"))
        self.Pic.setFrameShape(QtWidgets.QFrame.Box)
        self.Pic.setObjectName("Pic")
        self.Helping = QtWidgets.QLabel(self.Tab4)
        self.Helping.setGeometry(QtCore.QRect(125, 490, 721, 51))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(22)
        font.setItalic(True)
        self.Helping.setFont(font)
        self.Helping.setAlignment(QtCore.Qt.AlignCenter)
        self.Helping.setObjectName("Helping")
        self.Mail = QtWidgets.QLabel(self.Tab4)
        self.Mail.setGeometry(QtCore.QRect(110, 545, 721, 51))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(22)
        font.setItalic(True)
        self.Mail.setFont(font)
        self.Mail.setAlignment(QtCore.Qt.AlignCenter)
        self.Mail.setObjectName("Mail")
        
        self.Tabs.addTab(self.Tab4, "")
        
        self.stackedWidget.addWidget(self.Page_1)
        self.Page_2 = QtWidgets.QWidget()
        self.Page_2.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.Page_2.setObjectName("Page_2")
        self.Progress_Bar = QtWidgets.QProgressBar(self.Page_2)
        self.Progress_Bar.setGeometry(QtCore.QRect(250, 250, 541, 41))
        self.Progress_Bar.setCursor(QtGui.QCursor(QtCore.Qt.BusyCursor))
        self.Progress_Bar.setProperty("value", 0)
        self.Progress_Bar.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.Progress_Bar.setObjectName("Progress_Bar")
        self.Show_Results = QtWidgets.QPushButton(self.Page_2)
        self.Show_Results.hide()
        self.Show_Results.setGeometry(QtCore.QRect(440, 350, 161, 51))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(20)
        font.setItalic(True)
        self.Show_Results.setFont(font)
        self.Show_Results.setMouseTracking(False)
        self.Show_Results.setTabletTracking(False)
        self.Show_Results.setFocusPolicy(QtCore.Qt.NoFocus)
        self.Show_Results.setStyleSheet("background-color: rgb(204, 204, 204);\n"
"border-color: rgb(27, 27, 27);")
        self.Show_Results.setObjectName("Show_Results")
        
        self.Show_Results.clicked.connect(self.page_change2)
        
        self.stackedWidget.addWidget(self.Page_2)
        self.Page_3 = QtWidgets.QWidget()
        self.Page_3.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.Page_3.setObjectName("Page_3")
        self.Table = QtWidgets.QTableWidget(self.Page_3)
        self.Table.setGeometry(QtCore.QRect(40, 20, 700, 671))
        self.Table.setStyleSheet("background-color: rgb(255, 255, 255);")
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Table.sizePolicy().hasHeightForWidth())
        self.Table.setSizePolicy(sizePolicy)
        self.Table.setMinimumSize(QtCore.QSize(941, 671))
        self.Table.setMaximumSize(QtCore.QSize(941, 16777215))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(16)
        self.Table.setFont(font)
        self.Table.setMouseTracking(True)
        self.Table.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.Table.setLineWidth(1)
        self.Table.setMidLineWidth(2)
        self.Table.setAlternatingRowColors(True)
        self.Table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.Table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.Table.setGridStyle(QtCore.Qt.SolidLine)
        self.Table.setCornerButtonEnabled(True)
        self.Table.setRowCount(1)
        self.Table.setObjectName("Table")
        self.Table.setColumnCount(2)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        font = QtGui.QFont()
        font.setFamily("Z003 [urw]")
        font.setPointSize(20)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        item.setFont(font)
        self.Table.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        font = QtGui.QFont()
        font.setFamily("Z003 [urw]")
        font.setPointSize(20)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        item.setFont(font)
        self.Table.setHorizontalHeaderItem(1, item)


        self.Table.horizontalHeader().setVisible(True)
        self.Table.horizontalHeader().setCascadingSectionResizes(False)
        self.Table.horizontalHeader().setDefaultSectionSize(307)
        self.Table.horizontalHeader().setSortIndicatorShown(False)
        self.Table.horizontalHeader().setStretchLastSection(False)
        self.Table.verticalHeader().setVisible(True)
        self.Table.verticalHeader().setCascadingSectionResizes(False)
        self.Table.verticalHeader().setSortIndicatorShown(False)
        self.Table.verticalHeader().setStretchLastSection(False)
        

        
        self.stackedWidget.addWidget(self.Page_3)
        self.save_file = QtWidgets.QPushButton(self.Page_3)
        self.save_file.setGeometry(QtCore.QRect(650, 480, 75, 41))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(20)
        font.setItalic(True)
        self.save_file.setFont(font)
        self.save_file.setMouseTracking(False)
        self.save_file.setTabletTracking(False)
        self.save_file.setFocusPolicy(QtCore.Qt.NoFocus)
        self.save_file.setStyleSheet("background-color: rgb(204, 204, 204);\n"
"border-color: rgb(27, 27, 27);")
        self.save_file.setObjectName("save_file")
        
        self.save_file.clicked.connect(self.file_save)
        
        
        self.stackedWidget.addWidget(self.Page_3)
        self.show_chimera = QtWidgets.QPushButton(self.Page_3)
        self.show_chimera.setGeometry(QtCore.QRect(735, 480,150, 41))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(20)
        font.setItalic(True)
        self.show_chimera.setFont(font)
        self.show_chimera.setMouseTracking(False)
        self.show_chimera.setTabletTracking(False)
        self.show_chimera.setFocusPolicy(QtCore.Qt.NoFocus)
        self.show_chimera.setStyleSheet("background-color: rgb(204, 204, 204);\n"
"border-color: rgb(27, 27, 27);")
        self.show_chimera.setObjectName("show_chimera")  
        self.show_chimera.clicked.connect(self.chimera_command_generator)
#        self.show_chimera.clicked.connect(self.chimera_show)
        
        
        self.stackedWidget.addWidget(self.Page_3)
        self.show_plots = QtWidgets.QPushButton(self.Page_3)
        self.show_plots.setGeometry(QtCore.QRect(900, 480, 75, 41))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(20)
        font.setItalic(True)
        self.show_plots.setFont(font)
        self.show_plots.setMouseTracking(False)
        self.show_plots.setTabletTracking(False)
        self.show_plots.setFocusPolicy(QtCore.Qt.NoFocus)
        self.show_plots.setStyleSheet("background-color: rgb(204, 204, 204);\n"
"border-color: rgb(27, 27, 27);")
        self.show_plots.setObjectName("show_plots")  
        self.show_plots.clicked.connect(self.plots)
        self.Pic = QtWidgets.QLabel(self.Page_3)
        self.Pic.setGeometry(QtCore.QRect(560, 75, 400, 300))
        self.Pic.setText("")
        self.Pic.setStyleSheet("background-color: rgb(255,248,220);\n"
"border-color: rgb(27, 27, 27);")
        
        V_Motif.setCentralWidget(self.centralwidget)

        self.retranslateUi(V_Motif)
        QtCore.QMetaObject.connectSlotsByName(V_Motif)

    def retranslateUi(self, V_Motif):
        _translate = QtCore.QCoreApplication.translate
        V_Motif.setWindowTitle(_translate("V_Motif", "V_Motif"))
        self.VMotif.setText(_translate("V_Motif", "V_Motif"))
        self.VMotif_2.setText(_translate("V_Motif", "Conserved Patterns Identification"))
        self.Tabs.setTabText(self.Tabs.indexOf(self.Tab1), _translate("V_Motif", "Home"))
        self.L_Protein_ID.setText(_translate("V_Motif", "Enter PDB ID:"))
        self.B_Execute.setText(_translate("V_Motif", "Execute"))
        self.L_Output_F.setText(_translate("V_Motif", "Output Folder"))
        self.Tabs.setTabText(self.Tabs.indexOf(self.Tab2), _translate("V_Motif", "Finding Motifs"))
        self.ProDom.setText(_translate("V_Motif", "InterPro"))
        self.MEME.setText(_translate("V_Motif", "MEME"))
        self.PRINTS.setText(_translate("V_Motif", "MOTIF"))
        self.Databases.setText(_translate("V_Motif", "Following are different motif databases,"))
        self.Databases_2.setText(_translate("V_Motif", "Check them for any query or analysis."))
        self.Tabs.setTabText(self.Tabs.indexOf(self.Tab3), _translate("V_Motif", "Other Databases Link"))
        self.Helping.setText(_translate("V_Motif", "---->    For any other query, read out \"Helping Document\" "))
        self.Mail.setText(_translate("V_Motif", "---->    Contact at:   \"amanabatool@bs.qau.edu.pk\""))
        self.Tabs.setTabText(self.Tabs.indexOf(self.Tab4), _translate("V_Motif", "Help"))
        self.Show_Results.setText(_translate("V_Motif", "Show Results"))
        item = self.Table.horizontalHeaderItem(0)
        item.setText(_translate("V_Motif", "Motifs_"))
        item = self.Table.horizontalHeaderItem(1)
        item.setText(_translate("V_Motif", "Color in Chimera"))
        self.save_file.setText(_translate("V_Motif", "Save"))
        self.show_chimera.setText(_translate("V_Motif", "Visualize"))
        self.show_plots.setText(_translate("V_Motif", "Plots"))
        

    def page_change(self):
        global b
        global c
        global x
        
        self.get_ID()
        self.connect()
        
        print ("x = ", x)
        
        if (b != '') and  (c != '') and (x != ''):
            self.stackedWidget.setCurrentIndex(1)
            self.startss()
                
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("One or more Fields are empty or no Internet Connection")
            msg.setInformativeText('Fill out those fields or check your Internet Connection')
            msg.setWindowTitle("Error")
            msg.exec_()
        
        
    def page_change2(self):
        self.stackedWidget.setCurrentIndex(2)
        
    def Pro_Dom(self):
        webbrowser.open("https://www.ebi.ac.uk/interpro/search/sequence/")
        
    def ME_ME(self):
        webbrowser.open("https://meme-suite.org/meme/tools/meme")
        
    def PRI_NTS(self):
        webbrowser.open("https://www.genome.jp/tools/motif/")
        
    def startss(self):          
        self.t1=threading.Thread(target=self.code, args=())
        self.t1.start()
        #self.t7=threading.Thread(target=self.conf_file, args=())
        #self.t7.start()
            
    def get_out(self):
        global dir_name
        global b
        dir_name = QFileDialog.getExistingDirectory(None, 'Select Output Folder')
        if dir_name:
            self.Enter_Output.setText(dir_name)
            b = dir_name

        filenames = os.listdir (dir_name)
        global dir_name1
        dir_name1 = "Result"
        
    def get_ID(self):
        global idss
        global idsss
        global c
        global idd
        idss = self.Enter_ID.toPlainText()
        if idss:
            idsss = idss.lower()
            self.Enter_ID.setText(idsss)
            c = idsss
            idd = c
        _translate = QtCore.QCoreApplication.translate
        item = self.Table.horizontalHeaderItem(0)
        item.setText(_translate("V_Motif", "Motifs of {}".format(idss)))
        self.save_file.setText(_translate("V_Motif", "Save"))
        
    def code(self):
        self.pdb_id = c
        self.out_dir = dir_name1
        self.id_score = idscore

        self._motifs = None
        
        if self._motifs == None: 
            self.struct_extract()
            self.pdb_seq_extract()
            
            self.Progress_Bar.setProperty("value", 50)
            self.main()
            self.Progress_Bar.setProperty("value", 100)
        print ("Done in" , time.time())
        print ("Process Done")
        
        self.Show_Results.show()
    def connect(self):
        global connect
        connect = 1
        global x
        conn = httplib.HTTPConnection("www.google.com", timeout=3)
        try:
            conn.request("HEAD", "/")
            conn.close()
            print ("Available")
            x = connect
            print ("True")
        except Exception as e:
            print(e)
            print ("False")
            x = ''
            
        print (x)
        
    def struct_extract(self): 
        pdb_inst = StructSeqRetrieve(self.pdb_id, self.out_dir)
        pdb_inst.struct_retrieve()
        pdb_inst.replace_ent2pdb()

    def pdb_seq_extract(self): 
        pdb_seq = extract_seq(f"{self.out_dir}/{self.pdb_id}.pdb", "pdb-seqres")
        write_seq(pdb_seq, f"{self.out_dir}/{self.pdb_id}.fasta")
     
    def seq_frag(self):
        text_file = open(f"{self.out_dir}/{self.pdb_id}.fasta", "r+")
        new_f = text_file.readlines()
        text_file.seek(0)
        for line in new_f:
            if ">" not in line:
                text_file.write(line)
        text_file.truncate()
        text_file2 = open(f"{self.out_dir}/{self.pdb_id}.fasta", "r+")
        data = text_file2.read().replace('\n', '')
        n = 9
        new_input = '\n'.join(re.findall('.{1,%i}' % n, data))
        #print(new_input)
        seq_file = open('Seq_input.txt', 'w')
        seq_file.write(new_input)
        seq_file.close()
    def apply_matrix(self):
     #function to apply matrix on motif_list   
        global final_array
        final_array = []
        global count_none 
        count_none = 0
        global color_text 
        color_text = " "
        file = 'Seq_input.txt'
        
    
        new_motif = open_file1(file)
#        print("new motif" , new_motif)
    
        x = len(new_motif)
        
        for k in range ( x):
            
            if new_motif[k] == 'AR' :
                final_array.append("AR")
        
            elif new_motif[k] == 'AN' :
                final_array.append("AN")
        
            elif new_motif[k] == 'AD' :
                final_array.append("AD")
    
            elif new_motif[k] == 'AC' :
                final_array.append("AC")
        
            elif new_motif[k] == 'AE' :
                final_array.append("AE")
            
            elif new_motif[k] == 'AQ' :
                final_array.append("AQ")
        
            elif new_motif[k] == 'AG' :
                final_array.append("AG")
    
            elif new_motif[k] == 'AH' :
                final_array.append("AH")
        
            elif new_motif[k] == 'AI' :
                final_array.append("AI")
        
            elif new_motif[k] == 'AL' :
                final_array.append("AL")
        
            elif new_motif[k] == 'AK' :
                final_array.append("AK")
        
            elif new_motif[k] == 'AM' :
                final_array.append("AM")
        
            elif new_motif[k] == 'AF' :
                final_array.append("AF")
        
            elif new_motif[k] == 'AP' :
                final_array.append("AP")
        
            elif new_motif[k] == 'AS' :
                final_array.append("AS")
        
            elif new_motif[k] == 'AT' :
                final_array.append("AT")
        
            elif new_motif[k] == 'AW' :
                final_array.append("AW")
        
            elif new_motif[k] == 'AY' :
                final_array.append("AY")
        
            elif new_motif[k] == 'AV' :
                final_array.append("AV")
        
            elif new_motif[k] == 'AA' :
                final_array.append("AA")
        
            elif new_motif[k] == 'RN' :
                final_array.append("RN")
        
            elif new_motif[k] == 'RD' :
                final_array.append("RD")
        
            elif new_motif[k] == 'RC' :
                final_array.append("RC")
        
            elif new_motif[k] == 'RE' :
                final_array.append("RE")
        
            elif new_motif[k] == 'RQ' :
                final_array.append("RQ")
        
            elif new_motif[k] == 'RG' :
                final_array.append("RG")
        
            elif new_motif[k] == 'RH' :
                final_array.append("RH")
        
            elif new_motif[k] == 'RI' :
                final_array.append("RI")
        
            elif new_motif[k] == 'RL' :
                final_array.append("RL")
        
            elif new_motif[k] == 'RK' :
                final_array.append("RK")
        
            elif new_motif[k] == 'RM' :
                final_array.append("RM")
        
            elif new_motif[k] == 'RF' :
                final_array.append("RF")
        
            elif new_motif[k] == 'RP' :
                final_array.append("RP")
        
            elif new_motif[k] == 'RS' :
                final_array.append("RS")
        
            elif new_motif[k] == 'RT' :
                final_array.append("RT")
        
            elif new_motif[k] == 'RW' :
                final_array.append("RW")
        
            elif new_motif[k] == 'RY' :
                final_array.append("RY")
        
            elif new_motif[k] == 'RV' :
                final_array.append("RV")
        
            elif new_motif[k] == 'RA' :
                final_array.append("RA")
        
            elif new_motif[k] == 'RR' :
                final_array.append("RR")
        
            elif new_motif[k] == 'ND' :
                final_array.append("ND")
        
            elif new_motif[k] == 'NC' :
                final_array.append("NC")
        
            elif new_motif[k] == 'NE' :
                final_array.append("NE")
        
            elif new_motif[k] == 'NQ' :
                final_array.append("NQ")
        
            elif new_motif[k] == 'NG' :
                final_array.append("NG")
        
            elif new_motif[k] == 'NH' :
                final_array.append("NH")
        
            elif new_motif[k] == 'NI' :
                final_array.append("NI")
        
            elif new_motif[k] == 'NL' :
                final_array.append("NL")
        
            elif new_motif[k] == 'NK' :
                final_array.append("NK")
        
            elif new_motif[k] == 'NM' :
                final_array.append("NM")
        
            elif new_motif[k] == 'NF' :
                final_array.append("NF")
        
            elif new_motif[k] == 'NP' :
                final_array.append("NP")
        
            elif new_motif[k] == 'NS' :
                final_array.append("NS")
        
            elif new_motif[k] == 'NT' :
                final_array.append("NT")
        
            elif new_motif[k] == 'NW' :
                final_array.append("NW")
        
            elif new_motif[k] == 'NY' :
                final_array.append("NY")
        
            elif new_motif[k] == 'NV' :
                final_array.append("NV")
        
            elif new_motif[k] == 'NA' :
                final_array.append("NA")
        
            elif new_motif[k] == 'NR' :
                final_array.append("NR")
        
            elif new_motif[k] == 'NN' :
                final_array.append("NN")
        
            elif new_motif[k] == 'DC' :
                final_array.append("DC")
            
            elif new_motif[k] == 'DE' :
                final_array.append("DE")
        
            elif new_motif[k] == 'DQ' :
                final_array.append("DQ")
        
            elif new_motif[k] == 'DG' :
                final_array.append("DG")
        
            elif new_motif[k] == 'DH' :
                final_array.append("DH")
        
            elif new_motif[k] == 'DI' :
                final_array.append("DI")
        
            elif new_motif[k] == 'DL' :
                final_array.append("DL")
        
            elif new_motif[k] == 'DK' :
                final_array.append("DK")
        
            elif new_motif[k] == 'DM' :
                final_array.append("DM")
        
            elif new_motif[k] == 'DF' :
                final_array.append("DF")
        
            elif new_motif[k] == 'DP' :
                final_array.append("DP")
        
            elif new_motif[k] == 'DS' :
                final_array.append("DS")
        
            elif new_motif[k] == 'DT' :
                final_array.append("DT")
        
            elif new_motif[k] == 'DW' :
                final_array.append("DW")
        
            elif new_motif[k] == 'DY' :
                final_array.append("DY")
        
            elif new_motif[k] == 'DV' :
                final_array.append("DV")
        
            elif new_motif[k] == 'DA' :
                final_array.append("DA")
        
            elif new_motif[k] == 'DR' :
                final_array.append("DR")
    
            elif new_motif[k] == 'DN' :
                final_array.append("DN")
        
            elif new_motif[k] == 'DD' :
                final_array.append("DD")
        
            elif new_motif[k] == 'CE' :
                final_array.append("CE")
        
            elif new_motif[k] == 'CQ' :
                final_array.append("CQ")
        
            elif new_motif[k] == 'CG' :
                final_array.append("CG")
        
            elif new_motif[k] == 'CH' :
                final_array.append("CH")
        
            elif new_motif[k] == 'CI' :
                final_array.append("CI")
        
            elif new_motif[k] == 'CL' :
                final_array.append("CL")
        
            elif new_motif[k] == 'CK' :
                final_array.append("CK")
        
            elif new_motif[k] == 'CM' :
                final_array.append("CM")
        
            elif new_motif[k] == 'CF' :
                final_array.append("CF")
        
            elif new_motif[k] == 'CP' :
                final_array.append("CP")
        
            elif new_motif[k] == 'CS' :
                final_array.append("CS")
        
            elif new_motif[k] == 'CT' :
                final_array.append("CT")
        
            elif new_motif[k] == 'CW' :
                final_array.append("CW")
        
            elif new_motif[k] == 'CY' :
                final_array.append("CY")
        
            elif new_motif[k] == 'CV' :
                final_array.append("CV")
        
            elif new_motif[k] == 'CA' :
                final_array.append("CA")
        
            elif new_motif[k] == 'CR' :
                final_array.append("CR")
        
            elif new_motif[k] == 'CN' :
                final_array.append("CN")
        
            elif new_motif[k] == 'CD' :
                final_array.append("CD")
        
            elif new_motif[k] == 'CC' :
                final_array.append("CC")
        
            elif new_motif[k] == 'EQ' :
                final_array.append("EQ")
        
            elif new_motif[k] == 'EG' :
                final_array.append("EG")
        
            elif new_motif[k] == 'EH' :
                final_array.append("EH")
        
            elif new_motif[k] == 'EI' :
                final_array.append("EI")
        
            elif new_motif[k] == 'EL' :
                final_array.append("EL")
        
            elif new_motif[k] == 'EK' :
                final_array.append("EK")
        
            elif new_motif[k] == 'EM' :
                final_array.append("EM")
        
            elif new_motif[k] == 'EF' :
                final_array.append("EF")
        
            elif new_motif[k] == 'EP' :
                final_array.append("EP")
        
            elif new_motif[k] == 'ES' :
                final_array.append("ES")
        
            elif new_motif[k] == 'ET' :
                final_array.append("ET")
        
            elif new_motif[k] == 'EW' :
                final_array.append("EW")
        
            elif new_motif[k] == 'EY' :
                final_array.append("EY")
        
            elif new_motif[k] == 'EV' :
                final_array.append("EV")
        
            elif new_motif[k] == 'EA' :
                final_array.append("EA")
        
            elif new_motif[k] == 'ER' :
                final_array.append("ER")
        
            elif new_motif[k] == 'EN' :
                final_array.append("EN")
        
            elif new_motif[k] == 'ED' :
                final_array.append("ED")
        
            elif new_motif[k] == 'EC' :
                final_array.append("EC")
        
            elif new_motif[k] == 'EE' :
                final_array.append("EE")
        
            elif new_motif[k] == 'QG' :
                final_array.append("QG")
        
            elif new_motif[k] == 'QH' :
                final_array.append("QH")
        
            elif new_motif[k] == 'QI' :
                final_array.append("QI")
        
            elif new_motif[k] == 'QL' :
                final_array.append("QL")
        
            elif new_motif[k] == 'QK' :
                final_array.append("QK")
        
            elif new_motif[k] == 'QM' :
                final_array.append("QM")
        
            elif new_motif[k] == 'QF' :
                final_array.append("QF")
        
            elif new_motif[k] == 'QP' :
                final_array.append("QP")
        
            elif new_motif[k] == 'QS' :
                final_array.append("QS")
        
            elif new_motif[k] == 'QT' :
                final_array.append("QT")
        
            elif new_motif[k] == 'QW' :
                final_array.append("QW")
        
            elif new_motif[k] == 'QY' :
                final_array.append("QY")
        
            elif new_motif[k] == 'QV' :
                final_array.append("QV")
        
            elif new_motif[k] == 'QA' :
                final_array.append("QA")
        
            elif new_motif[k] == 'QR' :
                final_array.append("QR")
        
            elif new_motif[k] == 'QN' :
                final_array.append("QN")
        
            elif new_motif[k] == 'QD' :
                final_array.append("QD")
        
            elif new_motif[k] == 'QC' :
                final_array.append("QC")
        
            elif new_motif[k] == 'QE' :
                final_array.append("QE")
        
            elif new_motif[k] == 'QQ' :
                final_array.append("QQ")
        
            elif new_motif[k] == 'GH' :
                final_array.append("GH")
        
            elif new_motif[k] == 'GI' :
                final_array.append("GI")
        
            elif new_motif[k] == 'GL' :
                final_array.append("GL")
        
            elif new_motif[k] == 'GK' :
                final_array.append("GK")
        
            elif new_motif[k] == 'GM' :
                final_array.append("GM")
        
            elif new_motif[k] == 'GF' :
                final_array.append("GF")
        
            elif new_motif[k] == 'GP' :
                final_array.append("GP")
        
            elif new_motif[k] == 'GS' :
                final_array.append("GS")
        
            elif new_motif[k] == 'GT' :
                final_array.append("GT")
        
            elif new_motif[k] == 'GW' :
                final_array.append("GW")
        
            elif new_motif[k] == 'GY' :
                final_array.append("GY")
        
            elif new_motif[k] == 'GV' :
                final_array.append("GV")
        
            elif new_motif[k] == 'GA' :
                final_array.append("GA")
        
            elif new_motif[k] == 'GR' :
                final_array.append("GR")
        
            elif new_motif[k] == 'GN' :
                final_array.append("GN")
        
            elif new_motif[k] == 'GD' :
                final_array.append("GD")
        
            elif new_motif[k] == 'GC' :
                final_array.append("GC")
            
            elif new_motif[k] == 'GE' :
                final_array.append("GE")
        
            elif new_motif[k] == 'GQ' :
                final_array.append("GQ")
        
            elif new_motif[k] == 'GG' :
                final_array.append("GG")
        
            elif new_motif[k] == 'HI' :
                final_array.append("HI")
        
            elif new_motif[k] == 'HL' :
                final_array.append("HL")
        
            elif new_motif[k] == 'HK' :
                final_array.append("HK")
        
            elif new_motif[k] == 'HM' :
                final_array.append("HM")
        
            elif new_motif[k] == 'HF' :
                final_array.append("HF")
        
            elif new_motif[k] == 'HP' :
                final_array.append("HP")
        
            elif new_motif[k] == 'HS' :
                final_array.append("HS")
        
            elif new_motif[k] == 'HT' :
                final_array.append("HT")
        
            elif new_motif[k] == 'HW' :
                final_array.append("HW")
        
            elif new_motif[k] == 'HY' :
                final_array.append("HY")
        
            elif new_motif[k] == 'HV' :
                final_array.append("HV")
        
            elif new_motif[k] == 'HA' :
                final_array.append("HA")
        
            elif new_motif[k] == 'HR' :
                final_array.append("HR")
        
            elif new_motif[k] == 'HN' :
                final_array.append("HN")
        
            elif new_motif[k] == 'HD' :
                final_array.append("HD")
        
            elif new_motif[k] == 'HC' :
                final_array.append("HC")
        
            elif new_motif[k] == 'HE' :
                final_array.append("HE")
        
            elif new_motif[k] == 'HQ' :
                final_array.append("HQ")
        
            elif new_motif[k] == 'HG' :
                final_array.append("HG")
        
            elif new_motif[k] == 'HH' :
                final_array.append("HH")
        
            elif new_motif[k] == 'IL' :
                final_array.append("IL")
        
            elif new_motif[k] == 'IK' :
                final_array.append("IK")
        
            elif new_motif[k] == 'IM' :
                final_array.append("IM")
        
            elif new_motif[k] == 'IF' :
                final_array.append("IF")
        
            elif new_motif[k] == 'IP' :
                final_array.append("IP")
        
            elif new_motif[k] == 'IS' :
                final_array.append("IS")
        
            elif new_motif[k] == 'IT' :
                final_array.append("IT")
        
            elif new_motif[k] == 'IW' :
                final_array.append("IW")
        
            elif new_motif[k] == 'IY' :
                final_array.append("IY")
        
            elif new_motif[k] == 'IV' :
                final_array.append("IV")
        
            elif new_motif[k] == 'IA' :
                final_array.append("IA")
        
            elif new_motif[k] == 'IR' :
                final_array.append("IR")
        
            elif new_motif[k] == 'IN' :
                final_array.append("IN")
        
            elif new_motif[k] == 'ID' :
                final_array.append("ID")
        
            elif new_motif[k] == 'IC' :
                final_array.append("IC")
        
            elif new_motif[k] == 'IE' :
                final_array.append("IE")
        
            elif new_motif[k] == 'IQ' :
                final_array.append("IQ")
        
            elif new_motif[k] == 'IG' :
                final_array.append("IG")
        
            elif new_motif[k] == 'IH' :
                final_array.append("IH")
        
            elif new_motif[k] == 'II' :
                final_array.append("II")
        
            elif new_motif[k] == 'LK' :
                final_array.append("LK")
        
            elif new_motif[k] == 'LM' :
                final_array.append("LM")
        
            elif new_motif[k] == 'LF' :
                final_array.append("LF")
        
            elif new_motif[k] == 'LP' :
                final_array.append("LP")
        
            elif new_motif[k] == 'LS' :
                final_array.append("LS")
        
            elif new_motif[k] == 'LT' :
                final_array.append("LT")
        
            elif new_motif[k] == 'LW' :
                final_array.append("LW")
        
            elif new_motif[k] == 'LY' :
                final_array.append("LY")
        
            elif new_motif[k] == 'LV' :
                final_array.append("LV")
        
            elif new_motif[k] == 'LA' :
                final_array.append("LA")
        
            elif new_motif[k] == 'LR' :
                final_array.append("LR")
        
            elif new_motif[k] == 'LN' :
                final_array.append("LN")
        
            elif new_motif[k] == 'LD' :
                final_array.append("LD")
        
            elif new_motif[k] == 'LC' :
                final_array.append("LC")
        
            elif new_motif[k] == 'LE' :
                final_array.append("LE")
        
            elif new_motif[k] == 'LQ' :
                final_array.append("LQ")
        
            elif new_motif[k] == 'LG' :
                final_array.append("LG")
        
            elif new_motif[k] == 'LH' :
                final_array.append("LH")
        
            elif new_motif[k] == 'LI' :
                final_array.append("LI")
        
            elif new_motif[k] == 'LL' :
                final_array.append("LL")
        
            elif new_motif[k] == 'KM' :
                final_array.append("KM")
        
            elif new_motif[k] == 'KF' :
                final_array.append("KF")
        
            elif new_motif[k] == 'KP' :
                final_array.append("KP")
        
            elif new_motif[k] == 'KS' :
                final_array.append("KS")
        
            elif new_motif[k] == 'KT' :
                final_array.append("KT")
        
            elif new_motif[k] == 'KW' :
                final_array.append("KW")
        
            elif new_motif[k] == 'KY' :
                final_array.append("KY")
        
            elif new_motif[k] == 'KV' :
                final_array.append("KV")
        
            elif new_motif[k] == 'KA' :
                final_array.append("KA")
        
            elif new_motif[k] == 'KR' :
                final_array.append("KR")
        
            elif new_motif[k] == 'KN' :
                final_array.append("KN")
        
            elif new_motif[k] == 'KD' :
                final_array.append("KD")
        
            elif new_motif[k] == 'KC' :
                final_array.append("KC")
        
            elif new_motif[k] == 'KE' :
                final_array.append("KE")
            
            elif new_motif[k] == 'KQ' :
                final_array.append("KQ")
        
            elif new_motif[k] == 'KG' :
                final_array.append("KG")
            
            elif new_motif[k] == 'KH' :
                final_array.append("KH")
        
            elif new_motif[k] == 'KI' :
                final_array.append("KI")
        
            elif new_motif[k] == 'KL' :
                final_array.append("KL")
        
            elif new_motif[k] == 'KK' :
               final_array.append("KK")
    
            elif new_motif[k] == 'MF' :
                final_array.append("MF")
        
            elif new_motif[k] == 'MP' :
                final_array.append("MP")
        
            elif new_motif[k] == 'MS' :
                final_array.append("MS")
        
            elif new_motif[k] == 'MT' :
                final_array.append("MT")
        
            elif new_motif[k] == 'MW' :
                final_array.append("MW")
        
            elif new_motif[k] == 'MY' :
                final_array.append("MY")
        
            elif new_motif[k] == 'MV' :
                final_array.append("MV")
        
            elif new_motif[k] == 'MA' :
                final_array.append("MA")
        
            elif new_motif[k] == 'MR' :
                final_array.append("MR")
        
            elif new_motif[k] == 'MN' :
                final_array.append("MN")
        
            elif new_motif[k] == 'MD' :
                final_array.append("MD")
        
            elif new_motif[k] == 'MC' :
                final_array.append("MC")
        
            elif new_motif[k] == 'ME' :
                final_array.append("ME")
        
            elif new_motif[k] == 'MQ' :
                final_array.append("MQ")
        
            elif new_motif[k] == 'MG' :
                final_array.append("MG")
        
            elif new_motif[k] == 'MH' :
                final_array.append("MH")
        
            elif new_motif[k] == 'MI' :
                final_array.append("MI")
        
            elif new_motif[k] == 'ML' :
                final_array.append("ML")
        
            elif new_motif[k] == 'MK' :
                final_array.append("MK")
        
            elif new_motif[k] == 'MM' :
                final_array.append("MM")
        
            elif new_motif[k] == 'FP' :
                final_array.append("FP")
        
            elif new_motif[k] == 'FS' :
                final_array.append("FS")
        
            elif new_motif[k] == 'FT' :
                final_array.append("FT")
        
            elif new_motif[k] == 'FW' :
                final_array.append("FW")
        
            elif new_motif[k] == 'FY' :
                final_array.append("FY")
        
            elif new_motif[k] == 'FV' :
                final_array.append("FV")
        
            elif new_motif[k] == 'FA' :
                final_array.append("FA")
        
            elif new_motif[k] == 'FR' :
                final_array.append("FR")
        
            elif new_motif[k] == 'FN' :
                final_array.append("FN")
        
            elif new_motif[k] == 'FD' :
                final_array.append("FD")
        
            elif new_motif[k] == 'FC' :
                final_array.append("FC")
        
            elif new_motif[k] == 'FE' :
                final_array.append("FE")
        
            elif new_motif[k] == 'FQ' :
                final_array.append("FQ")
        
            elif new_motif[k] == 'FG' :
                final_array.append("FG")
        
            elif new_motif[k] == 'FH' :
                final_array.append("FH")
        
            elif new_motif[k] == 'FI' :
                final_array.append("FI")
        
            elif new_motif[k] == 'FL' :
                final_array.append("FL")
        
            elif new_motif[k] == 'FK' :
                final_array.append("FK")
        
            elif new_motif[k] == 'FM' :
                final_array.append("FM")
        
            elif new_motif[k] == 'FF' :
                final_array.append("FF")
        
            elif new_motif[k] == 'PS' :
                final_array.append("PS")
        
            elif new_motif[k] == 'PT' :
                final_array.append("PT")
        
            elif new_motif[k] == 'PW' :
                final_array.append("PW")
        
            elif new_motif[k] == 'PY' :
                final_array.append("PY")
            
            elif new_motif[k] == 'PV' :
                final_array.append("PV")
        
            elif new_motif[k] == 'PA' :
                final_array.append("PA")
        
            elif new_motif[k] == 'PR' :
                final_array.append("PR")
        
            elif new_motif[k] == 'PN' :
                final_array.append("PN")
        
            elif new_motif[k] == 'PD' :
                final_array.append("PD")
        
            elif new_motif[k] == 'PC' :
                final_array.append("PC")
        
            elif new_motif[k] == 'PE' :
                final_array.append("PE")
        
            elif new_motif[k] == 'PQ' :
                final_array.append("PQ")
        
            elif new_motif[k] == 'PG' :
                final_array.append("PG")
        
            elif new_motif[k] == 'PH' :
                final_array.append("PH")
        
            elif new_motif[k] == 'PI' :
                final_array.append("PI")
        
            elif new_motif[k] == 'PL' :
                final_array.append("PL")
        
            elif new_motif[k] == 'PK' :
                final_array.append("PK")
        
            elif new_motif[k] == 'PM' :
                final_array.append("PM")
        
            elif new_motif[k] == 'PF' :
                final_array.append("PF")
        
            elif new_motif[k] == 'PP' :
                final_array.append("PP")
        
            elif new_motif[k] == 'ST' :
                final_array.append("ST")
        
            elif new_motif[k] == 'SW' :
                final_array.append("SW")
        
            elif new_motif[k] == 'SY' :
                final_array.append("SY")
        
            elif new_motif[k] == 'SV' :
                final_array.append("SV")
        
            elif new_motif[k] == 'SA' :
                final_array.append("SA")
        
            elif new_motif[k] == 'SR' :
                final_array.append("SR")
        
            elif new_motif[k] == 'SN' :
                final_array.append("SN")
        
            elif new_motif[k] == 'SD' :
                final_array.append("SD")
        
            elif new_motif[k] == 'SC' :
                final_array.append("SC")
        
            elif new_motif[k] == 'SE' :
                final_array.append("SE")
        
            elif new_motif[k] == 'SQ' :
                final_array.append("SQ")
        
            elif new_motif[k] == 'SG' :
                final_array.append("SG")
        
            elif new_motif[k] == 'SH' :
                final_array.append("SH")
        
            elif new_motif[k] == 'SI' :
                final_array.append("SI")
        
            elif new_motif[k] == 'SL' :
                final_array.append("SL")
        
            elif new_motif[k] == 'SK' :
                final_array.append("SK")
        
            elif new_motif[k] == 'SM' :
                final_array.append("SM")
        
            elif new_motif[k] == 'SF' :
                final_array.append("SF")
        
            elif new_motif[k] == 'SP' :
                final_array.append("SP")
        
            elif new_motif[k] == 'SS' :
                final_array.append("SS")
        
            elif new_motif[k] == 'TW' :
                final_array.append("TW")
        
            elif new_motif[k] == 'TY' :
                final_array.append("TY")
        
            elif new_motif[k] == 'TV' :
                final_array.append("TV")
        
            elif new_motif[k] == 'TA' :
                final_array.append("TA")
        
            elif new_motif[k] == 'TR' :
                final_array.append("TR")
        
            elif new_motif[k] == 'TN' :
                final_array.append("TN")
        
            elif new_motif[k] == 'TD' :
                final_array.append("TD")
        
            elif new_motif[k] == 'TC' :
                final_array.append("TC")
        
            elif new_motif[k] == 'TE' :
                final_array.append("TE")
        
            elif new_motif[k] == 'TQ' :
                final_array.append("TQ")
        
            elif new_motif[k] == 'TG' :
                final_array.append("TG")
                
            elif new_motif[k] == 'TH' :
                final_array.append("TH")
        
            elif new_motif[k] == 'TI' :
                final_array.append("TI")
        
            elif new_motif[k] == 'TL' :
                final_array.append("TL")
        
            elif new_motif[k] == 'TK' :
                final_array.append("TK")
        
            elif new_motif[k] == 'TM' :
                final_array.append("TM")
        
            elif new_motif[k] == 'TF' :
                final_array.append("TF")
        
            elif new_motif[k] == 'TP' :
                final_array.append("TP")
        
            elif new_motif[k] == 'TS' :
                final_array.append("TS")
                
            elif new_motif[k] == 'TT' :
                final_array.append("TT")
        
            elif new_motif[k] == 'WY' :
                final_array.append("WY")
        
            elif new_motif[k] == 'WV' :
                final_array.append("WV")
        
            elif new_motif[k] == 'WA' :
                final_array.append("WA")
        
            elif new_motif[k] == 'WR' :
                final_array.append("WR")
        
            elif new_motif[k] == 'WN' :
                final_array.append("WN")
        
            elif new_motif[k] == 'WD' :
                final_array.append("WD")
        
            elif new_motif[k] == 'WC' :
                final_array.append("WC")
        
            elif new_motif[k] == 'WE' :
                final_array.append("WE")
        
            elif new_motif[k] == 'WQ' :
                final_array.append("WQ")
        
            elif new_motif[k] == 'WG' :
                final_array.append("WG")
        
            elif new_motif[k] == 'WH' :
                final_array.append("WH")
        
            elif new_motif[k] == 'WI' :
                final_array.append("WI")
        
            elif new_motif[k] == 'WL' :
                final_array.append("WL")
        
            elif new_motif[k] == 'WK' :
                final_array.append("WK")
        
            elif new_motif[k] == 'WM' :
                final_array.append("WM")
        
            elif new_motif[k] == 'WF' :
                final_array.append("WF")
        
            elif new_motif[k] == 'WP' :
                final_array.append("WP")
        
            elif new_motif[k] == 'WS' :
                final_array.append("WS")
        
            elif new_motif[k] == 'WT' :
                final_array.append("WT")
        
            elif new_motif[k] == 'WW' :
                final_array.append("WW")
        
            elif new_motif[k] == 'YV' :
                final_array.append("YV")
        
            elif new_motif[k] == 'YA' :
                final_array.append("YA")
        
            elif new_motif[k] == 'YR' :
                final_array.append("YR")
        
            elif new_motif[k] == 'YN' :
                final_array.append("YN")
        
            elif new_motif[k] == 'YD' :
                final_array.append("YD")
        
            elif new_motif[k] == 'YC' :
                final_array.append("YC")
        
            elif new_motif[k] == 'YE' :
                final_array.append("YE")
        
            elif new_motif[k] == 'YQ' :
                final_array.append("YQ")
        
            elif new_motif[k] == 'YG' :
                final_array.append("YG")
        
            elif new_motif[k] == 'YH' :
                final_array.append("YH")
                
            elif new_motif[k] == 'YI' :
                final_array.append("YI")
        
            elif new_motif[k] == 'YL' :
                final_array.append("YL")
        
            elif new_motif[k] == 'YK' :
                final_array.append("YK")
        
            elif new_motif[k] == 'YM' :
                final_array.append("YM")
        
            elif new_motif[k] == 'YF' :
                final_array.append("YF")
        
            elif new_motif[k] == 'YP' :
                final_array.append("YP")
        
            elif new_motif[k] == 'YS' :
                final_array.append("YS")
        
            elif new_motif[k] == 'YT' :
                final_array.append("YT")
        
            elif new_motif[k] == 'YW' :
                final_array.append("YW")
                
            elif new_motif[k] == 'YY' :
                final_array.append("YY")
        
            elif new_motif[k] == 'VA' :
                final_array.append("VA")
        
            elif new_motif[k] == 'VR' :
                final_array.append("VR")
            
            elif new_motif[k] == 'VN' :
                final_array.append("VN")
        
            elif new_motif[k] == 'VD' :
                final_array.append("VD")
        
            elif new_motif[k] == 'VC' :
                final_array.append("VC")
        
            elif new_motif[k] == 'VE' :
                final_array.append("VE")
        
            elif new_motif[k] == 'VQ' :
                final_array.append("VQ")
        
            elif new_motif[k] == 'VG' :
                final_array.append("VG")
        
            elif new_motif[k] == 'VH' :
                final_array.append("VH")
        
            elif new_motif[k] == 'VI' :
                final_array.append("VI")
        
            elif new_motif[k] == 'VL' :
                final_array.append("VL")
        
            elif new_motif[k] == 'VK' :
                final_array.append("VK")
        
            elif new_motif[k] == 'VM' :
                final_array.append("VM")
        
            elif new_motif[k] == 'VF' :
                final_array.append("VF")
        
            elif new_motif[k] == 'VP' :
                final_array.append("VP")
        
            elif new_motif[k] == 'VS' :
                final_array.append("VS")
        
            elif new_motif[k] == 'VT' :
                final_array.append("VT")
        
            elif new_motif[k] == 'VW' :
                final_array.append("VW")
        
            elif new_motif[k] == 'VY' :
                final_array.append("VY")
        
            elif new_motif[k] == 'VV' :
                final_array.append("VV")
    
            else:
                count_none = count_none+1
    
    
            k = k+1       
        
        
        print (final_array)
        df_motif = pd.DataFrame(final_array)
        df_motif.to_csv('AA_file_without_freq.csv')
           
     
        
    def main(self):
        self.seq_frag()
        self.apply_matrix()
        df = pd.read_csv('AA_file_without_freq.csv', header=None)
        #open the csv file generated by V_motif_database.py
        df.rename(columns={0: 'Sr. No.', 1: 'Neighbouring_residues'}, inplace=True)
        #change the name of columns
        df.to_csv('myfile_with_col.csv', index=False)
        neighbor = "Neighbouring_residues"
        sr_no = 'Sr. No.'
        df_1 = pd.read_csv('freq.csv')
        df_2 = pd.read_csv('myfile_with_col.csv')
    
        df_1[neighbor] = df_1[neighbor].str.strip()
    
        df_3 = pd.merge(df_1, df_2, on=neighbor)
    
        df_3.drop(sr_no, axis=1, inplace=True)
    
        df_3.to_csv('merged.csv', index=False)
        df_final = pd.read_csv('merged.csv')
        print(df_final.loc[df_final['Frequency'] >= 0.6].to_csv('threshold_frequency.csv',index=False))
        print("Residues saved in file:threshold_frequency")
        print(df_final)
        seq_file2 = open('Seq_input.txt','r')
        seq_1 = seq_file2.read()
        text_file3 = open("data_seq.txt", "w") 
        for i in seq_1:
            seq_1 = seq_1[1:] +  "\n"
        #write string to file
            text_file3.write(seq_1)
        file = "data_seq.txt"
        freq_file = 'threshold_frequency.csv'
        resid_col = 'Neighbouring_residues'
        out_file = 'predicted_motifs_from_matrix_application.txt'
        
            
        data = self.read_file(file)
        df_new = pd.read_csv(freq_file)
        patterns = df_new[resid_col].tolist()
        found = self.find_pattern(data, patterns)
        self.write_to_file(found, out_file)
        self.Table.resizeRowsToContents()
        with open('predicted_motifs_from_matrix_application.txt', 'r') as f:
            lines = f.read().splitlines()
    
        lengths = (len(lines))
        self.Table.setRowCount(lengths)
        rgb_colors = ["green","orangered","violet","red","teal","deeppink","yellow","brown","blue","cyan","gray","plum","khaki","salmon","fuchsia" ,"turquoise","lightgreen","indigo","goldenrod","thistle","paleturquoise","antiquewhite","forestgreen","olivedrab","tomato","rosybrown","skyblue","mediumpurple","mediumvioletred","azure","aquamarine","darkmagenta","darkslateblue","deepskyblue","lightseagreen","peru ","coral","firebrick","orchid","papayawhip","limegreen","darkorange","sandybrown","mistyrose","darkorchid","sienna ","darkkhaki","lavendergreen"]
        j = 0
        new_len = (len(rgb_colors))
        i = 0
        while i < lengths:            
                item = QtWidgets.QTableWidgetItem(lines[i])
                item.setBackground(QtGui.QColor(255,255,255))
                self.Table.setItem(i, 0, item)
                i += 1
                
        while j < new_len:            
                item_1 = QtWidgets.QTableWidgetItem(rgb_colors[j])
                item_1.setBackground(QtGui.QColor(255,255,255))
                self.Table.setItem(j, 1, item_1)
                j += 1
            
                    
                    
    def file_save(self):
            path, ok = QtWidgets.QFileDialog.getSaveFileName(
                None, 'Save CSV', os.getenv('HOME'), 'CSV(*.csv)')
            if ok:
                columns = range(self.Table.columnCount())
                header = [self.Table.horizontalHeaderItem(column).text()
                          for column in columns]
                with open(path, 'w') as csvfile:
                    writer = csv.writer(
                        csvfile, dialect='excel', lineterminator='\n')
                    writer.writerow(header)
                    for row in range(self.Table.rowCount()):
                        writer.writerow(
                            self.Table.item(row, column).text()
                            for column in columns)
            print("File saved")                                
        
    def write_to_file(self,data, dir):
        with open(dir, "w") as file:
            for line in data:
                print(line, file=file)
    
    def read_file(self,dir):
        with open(dir, "r") as file:
            data = file.read().split()
    
        data = [pat for pat in data if len(pat) == 9]
        return data
   
    def find_pattern(self,total_list, pattern_list):
        found = []
        for motif in total_list:
            for three_comb in itertools.combinations(pattern_list, 3):
                first_comb, second_comb, third_comb = three_comb
                if first_comb in motif and second_comb in motif and third_comb in motif:
                    found.append(motif)
        return set(found) 

    def open_file(self, file_name):
        dir_name=f"{self.out_dir}/{self.pdb_id}.pdb"
        traj=mdt.load(dir_name)
        top=traj.topology
        df, bounds = top.to_dataframe()
        df.to_csv("Orinnal.csv")
        return df
    
    def seq_generator(self):
        df=self.open_file(f"{self.out_dir}/{self.pdb_id}.pdb")
        sequence=[]
        chains=df['chainID'].unique()
        for c in chains:
            names=df[df['chainID']==c]
            names=names.drop_duplicates(subset="resSeq")
            n=names['resName']
            n=n.to_list()
            n=str(n)
            n=n.replace('[',"")
            n=n.replace(']',"")
            n=n.replace(',',"")
            n=n.replace("'","")
            n=n.replace(' ',"")
            n=seq1(n)
            sequence.append(str(n))
        return sequence
    def chimera_command_generator(self):
        seq=self.seq_generator()
        file = open("predicted_motifs_from_matrix_application.txt")
        motifs=file.readlines()
        mff=[]
        for m in motifs:
            a=m.replace("\n","")
            mff.append(a)
        chains = 'a b c d e f g h i j k l m n o p q r s t u v w x y z'.split()
        file=open('chimera_colors.txt')
        colors=file.readlines()
        file.close()
        chimera_file=open("chimera_script.py","w+")
        chimera_file.write("from chimera import openModels\nfrom chimera import runCommand as rc")
        chimera_file.write(f"\nopenModels.open('{self.out_dir}/{self.pdb_id}.pdb')")
        for s in range(len(seq)):
            count =0
            for m in mff:
                count = count+1
                r=re.search(m, seq[s])
                
                if r != None:
                    temp=np.arange(r.span()[0]+1, r.span()[1]+1)
                    temp = list(temp)
                    if len(temp) > 0:
                        chimera_color=str(f"color {colors[count]} : ".replace("\n","") + str(temp).replace('[',"").replace(']',"").replace(",", f".{chains[s]}, ")+f".{chains[s]}")
                        a=f"color {colors[count]} : ".replace("\n","") + str(temp).replace('[',"").replace(']',"")#.replace(" ", f".{chains[s]}, ")#+f".{chains[s]}"
                        chimera_file.write(f"\nrc('''{chimera_color}''')")
        
        chimera_file.close()
        os.system("/home/student/.local/UCSF-Chimera64-1.16/bin/chimera -- %F chimera_script.py")
    def plots(self):
        text_file_motifs = open("predicted_motifs_from_matrix_application.txt", "r")
        data_motif = text_file_motifs.read()
        text_file_motifs.close()
        List_of_hydrophobic_residues = ['G','A','V','L','I','P','F','M','W']
        List_of_polar_residues = ['S','T','C','N','Q','Y']
        List_of_Acidic_residues = ['D','E']
        List_of_Basic_residues = ['R','K','H'] 
        #print(data)
        #total_char = len(data_motif)
        print("PLOTS DISPLAYED")
        count_A = 0
        for i in data_motif:
            if i == "A":
                count_A = count_A + 1
        #print(count_A)
        count_G = 0
        for i in data_motif:
            if i == "G":
                count_G = count_G + 1
        #print(count_G)
        count_V = 0
        for i in data_motif:
            if i == "V":
                count_V = count_V + 1
        #print(count_V)
        count_L = 0
        for i in data_motif:
            if i == "L":
                count_L = count_L + 1
        #print(count_L)
        count_I = 0
        for i in data_motif:
            if i == "I":
                count_I = count_I + 1
        #print(count_I)
        count_P = 0
        for i in data_motif:
            if i == "P":
                count_P = count_P + 1
        #print(count_P)
        count_F = 0
        for i in data_motif:
            if i == "F":
                count_F = count_F + 1
        #print(count_F)
        count_M = 0
        for i in data_motif:
            if i == "M":
                count_M = count_M + 1
        #print(count_M)
        count_W = 0
        for i in data_motif:
            if i == "W":
                count_W = count_W + 1
        #print(count_W)
        count_S = 0
        for i in data_motif:
            if i == "S":
                count_S = count_S + 1
        #print(count_S)
        count_T = 0
        for i in data_motif:
            if i == "T":
                count_T = count_T + 1
        #print(count_T)
        count_C = 0
        for i in data_motif:
            if i == "C":
                count_C = count_C + 1
        #print(count_C)
        count_N = 0
        for i in data_motif:
            if i == "N":
                count_N = count_N + 1
        #print(count_N)
        count_Q = 0
        for i in data_motif:
            if i == "Q":
                count_Q = count_Q + 1
        #print(count_Q)
        count_Y = 0
        for i in data_motif:
            if i == "Y":
                count_Y = count_Y + 1
        #print(count_Y)
        count_D = 0
        for i in data_motif:
            if i == "D":
                count_D = count_D + 1
        #print(count_D)
        count_E = 0
        for i in data_motif:
            if i == "E":
                count_E = count_E + 1
        #print(count_E)
        count_R = 0
        for i in data_motif:
            if i == "R":
                count_R = count_R + 1
        #print(count_R)
        count_K = 0
        for i in data_motif:
            if i == "K":
                count_K = count_K + 1
        #print(count_K)
        count_H = 0
        for i in data_motif:
            if i == "H":
                count_H = count_H + 1
        #print(count_H)
        Hydrophobic_residues = count_A + count_G + count_V + count_L + count_I + count_P + count_F + count_M + count_W
        #print(Hydrophobic_residues)
        Polar_residues = count_S + count_T + count_C + count_N + count_Q + count_Y
        #print(Polar_residues)
        Acidic_residues = count_D + count_E
        #print(Acidic_residues)
        Basic_residues = count_R + count_K + count_H
        #print(Basic_residues)
        
        
        
        activities = ['Hydrophobic', 'Polar', 'Acidic', 'Basic']
     
        # portion covered by each label
        slices = [Hydrophobic_residues, Polar_residues, Acidic_residues, Basic_residues]
         
        # color for each label
        colors = ['r', 'y', 'g', 'b']
        plt.subplot(121) 
        # plotting the pie chart
        pie = plt.pie(slices, labels = activities, colors=colors,
                startangle=190, shadow = True, explode = (0, 0, 0.0, 0),
                radius = 1.0, autopct = '%1.1f%%')
         
        # plotting legend
        plt.legend(loc ="lower left", fontsize=5)
        left = [1, 2, 3, 4]
         
        # heights of bars
        height = [Hydrophobic_residues, Polar_residues, Acidic_residues, Basic_residues]
         
        # labels for bars
        tick_label = ['HP', 'Polar', 'Acidic', 'Basic']
        plt.subplot(122) 
        # plotting a bar chart
        bars = plt.bar(left, height, tick_label = tick_label,
                width = 0.7, color = ['red', 'yellow', 'green', 'blue'])
        plt.bar_label(bars, padding=0, fmt='%g')
         
        # naming the x-axis
        plt.xlabel('Properties', labelpad=7)
        # naming the y-axis
        plt.ylabel('No. of residues')
        # plot title
        plt.title('Physiochemical behaviour') 
        figure = plt.gcf()  # get current figure
        plt.savefig('filename.png')
        # function to show the plot
        plt.show()
        
        self.Pic.setPixmap(QtGui.QPixmap("filename.png"))
        self.Pic.setFrameShape(QtWidgets.QFrame.Box)
        self.Pic.setObjectName("Pic")
        





if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    V_Motif = QtWidgets.QMainWindow()
    ui = Ui_V_Motif()
    ui.setupUi(V_Motif)
    V_Motif.show()
    sys.exit(app.exec_())
