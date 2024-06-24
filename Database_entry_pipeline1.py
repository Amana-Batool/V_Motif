 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 12:26:50 2022

@author: student
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 13:32:59 2020

@author: nida
"""

import os
import glob
import shutil
import re

tool_files = os.getcwd()

import time
import threading

from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QTableWidgetItem
from PyQt5.QtWidgets import QFileDialog ,QMessageBox

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

from seq_retrieve import StructSeqRetrieve

from Bio.Blast import NCBIWWW

try:
    import httplib
except:
    import http.client as httplib

from aligner import xml_parser

from Bio.Align.Applications import ClustalOmegaCommandline

from read_write_seqs import extract_seq
from read_write_seqs import write_seq

from Stats import Analysis
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np

from identical_sequence_parser import IdenticalSequencesParser
from emboss import emboss_water
import subprocess

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
        self.VMotif.setGeometry(QtCore.QRect(340, 180, 261, 91))
        font = QtGui.QFont()
        font.setFamily("Z003 [urw]")
        font.setPointSize(50)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.VMotif.setFont(font)
        self.VMotif.setStyleSheet("background-color: rgb(102, 138, 139);")
        self.VMotif.setAlignment(QtCore.Qt.AlignCenter)
        self.VMotif.setObjectName("VMotif")
        self.VMotif_2 = QtWidgets.QLabel(self.Tab1)
        self.VMotif_2.setGeometry(QtCore.QRect(110, 290, 731, 61))
        font = QtGui.QFont()
        font.setFamily("Z003 [urw]")
        font.setPointSize(40)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
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
        self.B_Execute.setGeometry(QtCore.QRect(410, 320, 121, 51))
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
        self.L_Output_F.setGeometry(QtCore.QRect(220, 235, 211, 41))
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
        self.ProDom.setGeometry(QtCore.QRect(160, 345, 121, 61))
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
        self.MEME.setGeometry(QtCore.QRect(660, 345, 111, 61))
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
        self.PRINTS.setGeometry(QtCore.QRect(420, 345, 121, 61))
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
        self.Databases.setGeometry(QtCore.QRect(210, 155, 521, 51))
        font = QtGui.QFont()
        font.setFamily("Yrsa Light")
        font.setPointSize(22)
        font.setItalic(True)
        self.Databases.setFont(font)
        self.Databases.setAlignment(QtCore.Qt.AlignCenter)
        self.Databases.setObjectName("Databases")
        self.Databases_2 = QtWidgets.QLabel(self.Tab3)
        self.Databases_2.setGeometry(QtCore.QRect(210, 230, 521, 51))
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
        self.Pic.setPixmap(QtGui.QPixmap("Capture.png"))
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
        self.Table.setGeometry(QtCore.QRect(40, 20, 941, 671))
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
        self.Table.setLineWidth(2)
        self.Table.setMidLineWidth(2)
        self.Table.setAlternatingRowColors(True)
        self.Table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.Table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.Table.setGridStyle(QtCore.Qt.SolidLine)
        self.Table.setCornerButtonEnabled(True)
        self.Table.setRowCount(1)
        self.Table.setObjectName("Table")
        self.Table.setColumnCount(3)
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
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        font = QtGui.QFont()
        font.setFamily("Z003 [urw]")
        font.setPointSize(20)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        item.setFont(font)
        self.Table.setHorizontalHeaderItem(2, item)

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
        V_Motif.setCentralWidget(self.centralwidget)

        self.retranslateUi(V_Motif)
        QtCore.QMetaObject.connectSlotsByName(V_Motif)

    def retranslateUi(self, V_Motif):
        _translate = QtCore.QCoreApplication.translate
        V_Motif.setWindowTitle(_translate("V_Motif", "V_Motif"))
        self.VMotif.setText(_translate("V_Motif", "V_Motif"))
        self.VMotif_2.setText(_translate("V_Motif", "Identification of Conserved Patterns"))
        self.Tabs.setTabText(self.Tabs.indexOf(self.Tab1), _translate("V_Motif", "Home"))
        self.L_Protein_ID.setText(_translate("V_Motif", "Enter PDB ID:"))
        self.B_Execute.setText(_translate("V_Motif", "Execute"))
        self.L_Output_F.setText(_translate("V_Motif", "Output Folder"))
        self.Tabs.setTabText(self.Tabs.indexOf(self.Tab2), _translate("V_Motif", "Finding Motifs"))
        self.ProDom.setText(_translate("V_Motif", "ProDom"))
        self.MEME.setText(_translate("V_Motif", "MEME"))
        self.PRINTS.setText(_translate("V_Motif", "PRINTS"))
        self.Databases.setText(_translate("V_Motif", "Following are different motif databases,"))
        self.Databases_2.setText(_translate("V_Motif", "Check them for any query or analysis."))
        self.Tabs.setTabText(self.Tabs.indexOf(self.Tab3), _translate("V_Motif", "Other Databases Link"))
        self.Helping.setText(_translate("V_Motif", "---->    For any other query, read out \"Helping Document\" "))
        self.Mail.setText(_translate("V_Motif", "---->    Contact at:   \"nida.ansari@bs.qau.edu.pk\""))
        self.Tabs.setTabText(self.Tabs.indexOf(self.Tab4), _translate("V_Motif", "Help"))
        self.Show_Results.setText(_translate("V_Motif", "Show Results"))
        item = self.Table.horizontalHeaderItem(0)
        item.setText(_translate("V_Motif", "Motifs"))
        item = self.Table.horizontalHeaderItem(1)
        item.setText(_translate("V_Motif", "Repeating Protein IDs"))
        item = self.Table.horizontalHeaderItem(2)
        item.setText(_translate("V_Motif", "No. of Repeats"))

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
        webbrowser.open("http://prodom.prabi.fr/prodom/current/html/form.php")
        
    def ME_ME(self):
        webbrowser.open("https://meme-suite.org/meme/tools/meme")
        
    def PRI_NTS(self):
        webbrowser.open("http://130.88.97.239/PRINTS/index.php")
        
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
        dir_name1 = "//home//student//Documents//V_Motif_Results//V_Motif_Results_1"
        
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
        
    def code(self):
        self.pdb_id = c
        self.out_dir = dir_name1
        self.id_score = idscore

        self._motifs = None
        
        if self._motifs == None: 
            self.struct_extract()
            self.pdb_seq_extract()
            
            self.Progress_Bar.setProperty("value", 20)
            
            self.psi_blast()
            
            self.Progress_Bar.setProperty("value", 65)
            
            psi_results = self.parse_psi_xml()
            write_seq(psi_results, f"{self.out_dir}/seqs.fasta")

            self.water()
            self.seq_trimmer()
            
            self.Progress_Bar.setProperty("value", 80)

            write_seq(self.result_record, f"{self.out_dir}/trimmed_seqs.fasta")

            self.clustalomega()
            
            self.Progress_Bar.setProperty("value", 90)
            
            self.motif_finder()
            
            self.Progress_Bar.setProperty("value", 100)
            

            #except IOError as e: 
                #print("Could not find any meaningful target sequences.\nMotif finding failed: ", e)
              
        print ("Done in" , time.time())
        print ("Process Done")
        
        self.Show_Results.show()
        
        self.myfunc()
        
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

    def psi_blast(self): 
        print('Doing the BLAST and retrieving the results...')
        
        global x
        
        self.connect()
        
        print ("x = ", x)
        
        if (x != ''):
            print ("BLAST Start")
            fastafile = open(f"{self.out_dir}/{self.pdb_id}.fasta")
            result_handle = NCBIWWW.qblast('blastp', 'pdb', fastafile.read())
            fastafile.close()
        
            self.Progress_Bar.setProperty("value", 50)
        
            print ("BLAST End")

            blast_file = open(f"{self.out_dir}/psi.xml", "w")
            blast_file.write(result_handle.read())
            blast_file.close()
            result_handle.close()
        
            print ("Done BLAST")
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("No Internet Connection")
            msg.setInformativeText('Check your Internet Connection and Restart the Tool')
            msg.setWindowTitle("Error")
            msg.exec_()
        
        #file = (f"{self.out_dir}/{self.pdb_id}.fasta")

        #query = SeqIO.read(file, format="fasta")
        #result_handle = NCBIWWW.qblast("blastp", "pdb", query.seq)

    def parse_psi_xml(self):
        
        return xml_parser(f"{self.out_dir}/psi.xml")

    def water(self): 
        
        emboss_water(f"{self.out_dir}/{self.pdb_id}.fasta", f"{self.out_dir}/seqs.fasta", f"{self.out_dir}/water.fasta")

    def seq_trimmer(self): 
        file = f"{self.out_dir}/water.fasta"
        needle_record = list(AlignIO.parse(file, "emboss"))
        self.result_record = []

        for rec in needle_record[1:]: 
            reference_seq = rec[0]
            seq_parser = IdenticalSequencesParser(reference_seq, rec[1])

            result = seq_parser.highly_identical_seqs()
            if result:
                self.result_record.append(result)

    def clustalomega(self): 
        pathh = os.getcwd()
        os.chdir(tool_files)
        shutil.copy("c", f"{self.out_dir}")
        os.chdir(f"{self.out_dir}")

        in_file = "trimmed_seqs.fasta"
        out_file = "aligned_seq.fasta"
        clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
        cmd_str = str(clustalomega_cline).split(" ")
        cmd_str[0] = "./c"

        subprocess.run(cmd_str, check="True")

        os.remove("c")
        
        os.chdir(pathh)

    def motif_finder(self): 
        path = f"{self.out_dir}/Aligned_seq.fasta" 
        path2 = f"{self.out_dir}/aligned_seq.fasta"

        with open(path2) as f_input, open(path, 'w') as f_output:
            block = []

            for line in f_input:
                if line.startswith('>'):
                    if block:
                        f_output.write(''.join(block) + '\n')
                        block = []
                    f_output.write(line)
                else:
                    block.append(line.strip())

            if block:
                f_output.write(''.join(block) + '\n')
        
        os.remove(f"{self.out_dir}/aligned_seq.fasta")
        
        seqs = extract_seq(f"{self.out_dir}/Aligned_seq.fasta", "fasta")
        seq = [[x for x in y] for y in seqs]

        c = Analysis(seq)
        self.c_ent = c.conservation_score()
        self.norm_data = c.normalize_data(self.c_ent)
        
        self.norm_data_len = np.arange(1, len(self.c_ent) + 1)
        
        sns.lineplot(self.norm_data_len, self.norm_data)
        plt.title("Conservation score per amino acid position")
        plt.xlabel("Amino acid position")
        plt.ylabel("Normalized conservation score")
        plt.savefig(f"{self.out_dir}/Conservation_Plot.png")
        plt.show()
        
        cons_data = {num : data for num, data in zip(self.norm_data_len, self.norm_data)}
        c.csv_writer(f"{self.out_dir}/Motif_Score.csv", cons_data)
        
        conserved_posits = c.find_motif(self.norm_data, self.norm_data_len, threshold = 0)
        
        print (conserved_posits)
        
        lengths = (len(conserved_posits))
       
        self.Table.setRowCount(lengths)
        
        motif = ""
        id_list = ""
        n = 0
        m = 0
        
        wb = Workbook()
        out = wb.add_sheet('Sheet 1')
        
        out.write(0, 0, 'pdb_id')
        out.write(0, 1, 'Motifs')
        out.write(0, 2, 'Target_IDs')
        out.write(0, 3, 'Repeats')

        for key, val in conserved_posits.items():
            print (key, val)
            for x in range(key-1, val):
                y = seq[0][x]
                motif += (y)
                
            print (motif)
            
            if not re.match(r'^[_\W]+$', motif):
                self.Table.resizeRowsToContents()
                #self.Table.resizeColumnsToContents()
            
                item = QTableWidgetItem(motif) # create the item
                item.setTextAlignment(Qt.AlignHCenter) # change the alignment
                font = QtGui.QFont()
                font.setFamily("Yrsa Light")
                font.setPointSize(14)
                font.setItalic(True)
                item.setFont(font)
                self.Table.setItem(m, 0, item)
        
                infile = open(f"{self.out_dir}/Aligned_seq.fasta",'r')
                pattern = re.compile(motif)
                for line in infile:
                    line = line.strip("\n")
                    if line.startswith('>'):
                        name = line
                    else:      
                        s = re.findall(pattern, line)
                        if (s != []):
                            #print ('%s:%s' %(name,s))
                            id_list += name
                            id_list += " "
                            n += 1
                
                out.write(m+1, 0, idss)
                out.write(m+1, 1, motif)
                out.write(m+1, 2, id_list)
                out.write(m+1, 3, n)
            
                item = QTableWidgetItem(str(n)) # create the item
                item.setTextAlignment(Qt.AlignHCenter) # change the alignment
                font = QtGui.QFont()
                font.setFamily("Yrsa Light")
                font.setPointSize(14)
                font.setItalic(True)
                item.setFont(font)
                self.Table.setItem(m, 2, item)
                #self.table.setItem(m, 3, QTableWidgetItem(str(n)))
            
                item = QTableWidgetItem(id_list) # create the item
                item.setTextAlignment(Qt.AlignHCenter) # change the alignment
                font = QtGui.QFont()
                font.setFamily("Yrsa Light")
                font.setPointSize(14)
                font.setItalic(True)
                item.setFont(font)
                self.Table.setItem(m, 1, item)
                #self.table.setItem(m, 2, QTableWidgetItem(id_list))
                motif = ""
                id_list = ""
                n = 0
                m += 1
                
            else:
                continue

        wb.save(f"{self.out_dir}/Motif_List.xls")

        c.pymol_script_writer(f"{self.out_dir}/{self.pdb_id}_pymol.pml", conserved_posits, f"{self.out_dir}/{self.pdb_id}.pdb")
        
        #c.motif_coordinates(f"{self.out_dir}/{self.pdb_id}_motif_coordinates.txt", f"{self.out_dir}/{self.pdb_id}.pdb", 
        #conserved_posits)
          
    def joins(self): 
        self.t4.join()
        self.t5.join()
        self.t6.join()
        self.t1.join()
        self.t7.join()
        self.t8.join()
    def myfunc(self):
        import pandas as pd
        import xlrd
        import openpyxl as xl
        import mysql.connector
        import xlrd
        conn=mysql.connector.connect(host="localhost",user="root",passwd="welcome123",database="V_Motif_Database")
        cur=conn.cursor()
        
        loc = ("//home//student//Documents//V_Motif_Results//V_Motif_Results_1//Motif_List.xls")
        workbook = pd.read_excel(loc)
        l=list()
        workbook.head()
        #print(workbook)
        index = workbook.index
        number_of_rows = len(index)
        #print(number_of_rows)
        a = xlrd.open_workbook(loc)
        sheet = a.sheet_by_index(0)
        sheet.cell_value(0,0)
        for i in range (1,number_of_rows):
            l.append(tuple(sheet.row_values(i)))
        q = "INSERT INTO Motif_Table (PDB_Id, Motif, Target_Ids, Repeats) VALUES (%s, %s, %s, %s)"
        cur.executemany(q,l)
        conn.commit()
        print(cur.rowcount, "record inserted.")
        sql_select_Query = "select * from Motif_Table"
        cursor = conn.cursor()
        cursor.execute(sql_select_Query)
            # get all records
        records = cursor.fetchall()
        print("Total number of rows in table: ", cursor.rowcount)
        conn.commit()
    #    print("\nPrinting each row")
        a_file = open("//home//student//Documents//V_Motif_Results//V_Motif_Results_1//motif.txt", "w")
        n = ''
        for row in records:
            n = row[1]
            a_file.write(n + "\n")
        a_file.close()
        print("motif.txt file generated")
        conn.close()
        

       
    

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    V_Motif = QtWidgets.QMainWindow()
    ui = Ui_V_Motif()
    ui.setupUi(V_Motif)
    V_Motif.show()
    sys.exit(app.exec_())
