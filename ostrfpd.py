#!/usr/bin/python  
# encoding: utf-8
''' (Revised version for submission)
 #######   ######  ######## ########  ######## ########  ########     
##     ## ##    ##    ##    ##     ## ##       ##     ## ##     ##    
##     ## ##          ##    ##     ## ##       ##     ## ##     ##    
##     ##  ######     ##    ########  ######   ########  ##     ##    
##     ##       ##    ##    ##   ##   ##       ##        ##     ##    
##     ## ##    ##    ##    ##    ##  ##       ##        ##     ##    
 #######   ######     ##    ##     ## ##       ##        ########   Version 0.01
-------------------------------------------------------------------
Filename : 'ostrfpd.py' [Single, all-in-one package script except for the primer3 plugins]
OSTRFPD  : Copyright (c) 2018 Mathema VB. All rights reserved.
Email    : vivek_mathema@hotmail.com

Citation: "Mathema VB, Dondorp AM, Imwong M (2019) 
           OSTRFPD: Multifunctional tool for genome-wide short tandem repeat analysis for DNA, 
           transcripts and amino acid sequences with integrated primer designer. Evolutionary 
           Bioinformatics.DOI: 10.1177/1176934319843130" 


(O)mni (S)hort (T)andem (R)epeat (F)inder & (P)rimer (D)esigner (OSTRFPD)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the  License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
-----------------------------------------------------------------------------------------
__author__      = "Faculty of Tropical Medicine, Mahidol Univeristy, Bangkok, thailand  |
__license__     = "GNU/GPL"                                                             |
__version__     = "0.01"                                                                |
__maintainer__  =  "VBM"                                                                |
__email__       = "(under review)"                                                      |
__status__      = "Production"                                                          |
-----------------------------------------------------------------------------------------

#---------------Important python Dependencies:
#  Python 3.0+  (v3.5.0 - v3.6.0 tested, Windows 10, Ubunutu 16.04 and Ubuntu 18.04.2 LTS tested)
#  biopython==1.72
#  PyQt5==5.9.1
#  future==0.16.0
#  Matplotlib==2.2.2 (optional) for drawing the barchart graphs of motif list (use specific version only if version compability issues appears)
#  Warnining! It is NOT RECOMMENDED to indlcue Matplot module if used for compiling and making standalone BINARIES of Ostrfpd for windows or LINUX
#  numpy==1.14.5     (optional) for drawing  the barchart (use specific version only if version compability issues appears)
#----------------
#  typical use (Command line interface):     python.exe ostrfpd.py -input FILENAME <-output FILENAME>
#  typical use (Graphical user interface):   python.exe ostrfpd.py -gui true 
#-------------------------------------------------------------

##TYPICAL STRUCTURE OF THE OSTRFPD FOLDER

OSTRFPD_v0.01     
|
|-----ostrfpd.py                     
|-----gui_interface.py
|-----ostrfpd.exe          [ Note: ‘ostrfpd’      for LINUX BINARIES supplied seperately ]
|-----primer3_core.exe     [ Note: ‘primer3_core’ for LINUX BINARIES supplied seperately ]
|-----run_ostrfpd_bin_gui.bat   
|-----run_ostrfpd_src_gui.bat
|-----dict_dna.txt
|-----dict_protein.txt
|-----input.fasta
|-----requirments.txt
|-----LICENCE
|-----rna_seq             (Note: subfolder)
     |-----ribosome1_seq_URS0000857690.fasta
     |-----ribosome2_seq_URS0000A6B25D.fasta

'''
#--------------------------------------------------------------
from __future__ import division
from __future__ import print_function  
from   PyQt5 import QtCore, QtGui, QtWidgets, uic
from   PyQt5.QtWidgets import QMessageBox, QFileDialog
from   gui_interface import Ui_MainWindow              
from   sys import version_info                                # for getting information of system 
from   datetime import datetime
from   itertools import combinations_with_replacement as cwr
from   Bio import SeqIO
from   Bio.Seq import Seq
from   Bio import pairwise2 
from   Bio.SeqUtils import GC
from   Bio.SeqUtils import seq3
import sys,subprocess
import gzip,os,sys                                            #  for reading .gz (gunzipped) files and normal os, input-output stuffs
import re, array                                              #  for  creating conventional arrays
import argparse, time                                         # for command lien argument passage and getting time-related infos 
import warnings                                               # for removing version rleated and otehr minor type-related warnings
#-------------------------------------------------------------
global fsc_uml                                                # global type vlaue  for recursive function
global fsc_umr                                                # global type vlaue  for recursive function
global tr_flag                                                # global Flag (boolean) type vlaue  for recursive function
global fsc_string                                             # MISA form of teh string for motif lenth and repeat number 
global dict_file                                              # name of dictionary file to output dictionary   
global primer_saveFileName                                    # to keep tarck of total primers made
global dict_warning                                           # to warn if dictionay file input has issues
#-------------------------------------------------------------

# GUI version module for basic error handeling and user frieldly INPUTS (PyQt5 related modules)
class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):            
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)          
        #-----------------------------------------------------Declare Buttons and connect                              
        self.Start_Scan.clicked.connect(self.GUI_module)
        self.PrimerFlag.clicked.connect(self.Activate_flanking_seq)
        self.rev_comp_flag.clicked.connect(self.check_rev_motif)        
        self.open_root.clicked.connect(self.open_folder)
        self.min_unit.valueChanged.connect(self.check_motif)
        self.max_unit.valueChanged.connect(self.check_motif)
        self.accurate_eng.clicked.connect(self.check_accurate_engine)
        self.default_eng.clicked.connect(self.check_default_engine)
        self.dict_mode.clicked.connect(self.check_dict_mode)
        self.dict_mode.clicked.connect(self.check_dict_range)
        self.Protein_button.clicked.connect(self.check_protein)
        self.min_fix.valueChanged.connect(self.check_repeat) 
        self.actionStartScan.triggered.connect(self.GUI_module)
        self.sourceFileDig.clicked.connect(self.set_sourcefilepath)
        self.resultFileDig.clicked.connect(self.set_resultfilepath)
        self.dictFileDig.clicked.connect(self.set_dictfilepath)
        self.primer_saveDlg.clicked.connect(self.set_primerfilepath)
        self.output_format.clicked.connect(self.check_rev_motif)
        self.about_button.clicked.connect(self.about_ostrfpd)
        #---------------------------------------------------------  

    def about_ostrfpd(self):
    	QMessageBox.about(self, "OSTRFPD v0.01", "Microsatellite search and primer designer!\n Developed by: Mathema VB, Dondorp AM, Imwong M (2019)")

    def check_rev_motif(self):
    	if (self.rev_comp_flag.isChecked()==True) and (self.output_format.currentItem().text().count("Default")== 0 or \
    		                                            self.Protein_button.isChecked()==True):
    		self.rev_comp_flag.setChecked(False)
    		QMessageBox.about(self, "OSTRFPD Warning!","Option not compatible for Protein sequence, Fasta or alignment output mode")

    def set_sourcefilepath(self):
    	fileName = self.openFileNameDialog() #if  self.openFileNameDialog() != None else ""
    	self.input_file.setText(fileName) if fileName != "" else None
    	self.output_file.setText(fileName +"_res.txt") if fileName  !="" else None
    	self.save_primer.setText(fileName +"_prm.txt") if fileName  !="" else None

    def set_resultfilepath (self):
    	fileName = self.saveFileDialog()
    	self.output_file.setText(fileName) if fileName  !="" else None

    def set_dictfilepath(self):
    	fileName = self.openFileNameDialog()
    	self.dict_string.setText(fileName) if fileName  != "" else None
        
    #----------------------------Ssave file dialogue box
    def set_primerfilepath(self):
        fileName = self.saveFileDialog()
        self.save_primer.setText(fileName) if fileName  !="" else None

    #-----------------------------Save file dialogue box 
    def openFileNameDialog(self):    
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Open sequence file DNA, RNA or Protein(.fasta, .fa, .gz) ","","All Files (*.*);;FASTA Files (*.fasta);Gunzip FASTAs (*.gz)", options=options)
        
        return fileName

    #-----------------------------Open file dialogue box
    def saveFileDialog(self):    
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Save OSTRFPD results","","All Files (*);;Text Files (*.txt)", options=options)

        return fileName

    #-----------------------------Checking the engine input
    def check_accurate_engine (self):
    	# only make the single tick for the engine
    	if (self.accurate_eng.isChecked() == True):
    		self.default_eng.setChecked(False)
    		self.dict_mode.setChecked(False)
    		self.check_motif()
    		return

    #----------------------------Checking the the default_engine input
    def check_default_engine (self):
    	# only make the single tick for the engine
    	if (self.default_eng.isChecked() == True):
    		self.accurate_eng.setChecked(False)
    		self.dict_mode.setChecked(False)
    		self.check_motif()
    		return

    # Checking the dictionary mode input
    def check_dict_mode (self):
    	# only make the single tick for the engine
    	if (self.dict_mode.isChecked() == True):
    		self.accurate_eng.setChecked(False)
    		self.default_eng.setChecked(False) 		
    		return

    def check_protein(self):
        
   		if (self.Protein_button.isChecked() == True):
   			self.rev_comp_flag.setChecked(False)           # remove teh reverse comeplent tag
   			self.check_motif()                             # check motif length

   		if self.PrimerFlag.isChecked() == True:
   			QMessageBox.about(self, "OSTRFPD Warning!", "The [Make primer] option is only valid for DNA sequence source")
   			self.PrimerFlag.setChecked(False)

   		if (self.min_fix.value() == 0):
   			self.misa_string.setText("5,4,3,3")
  
    def check_repeat(self):
   		# no reverse complements and motif length checks
   		if (self.min_fix.value()) == 0 and (self.misa_string.toPlainText() == ""):
   			QMessageBox.about(self, "Warning!", "Dictionay mode is only valid in Fixed Minimum Repeat mode.\
                                     MISA-formatted settings for repeat number should be set empty (blank) for using dictionary-based scan")
   			self.min_fix.setValue(15)

    def check_dict_range(self):
        if ((self.min_fix.value()) == 0 and (self.misa_string.toPlainText() != "")):

            response = QMessageBox.question(self,'Warning!', "Dictionay mode is only valid in Fixed Minimum Repeats.\
            	                                \nMISA-formatted settings for minimum repeat number MUST be set to empty (blank) while using this mode!.\
            	                                \n\nAre you sure to use Dictionary mode?.\
            	                                \nIf YES, Please change default values for unit motif length. Set the minimum repeat number accodring to your need.", QMessageBox.Yes | QMessageBox.No)
            
            tmp_min_unit =  int(self.min_unit.value())
            tmp_max_unit =  int(self.max_unit.value())
            tmp_misa     =  str(self.misa_string.toPlainText().strip())
            tmp_min_fix  =  int(self.min_fix.value())

            if response == QMessageBox.Yes:
                self.misa_string.setText("")
                self.min_fix.setValue (4)
                self.min_unit.setValue(8)
                self.max_unit.setValue(8)
                self.dict_mode.setChecked(True)
            else:
            	self.min_unit.setValue(tmp_min_unit)
            	self.max_unit.setValue(tmp_max_unit)
            	self.min_fix.setValue(tmp_min_fix)
            	self.misa_string.setText(tmp_misa)
            	self.dict_mode.setChecked(False)
            	self.accurate_eng.setChecked(False)
            	self.default_eng.setChecked(True)
            
            self.check_motif()
            return
    # Module for check the inputs in GUI mode        
    def check_motif(self):
    	# check if DNA or RNA sequence for Accurate engine selection
    	if ((self.DNA_button.isChecked() == True  or self.RNA_button.isChecked() == True ) and self.accurate_eng.isChecked() == True):
    		self.max_unit.setValue(6) if int(self.max_unit.value()) > 6 else int(self.max_unit.value())
    		self.min_unit.setValue(self.max_unit.value()) if int(self.min_unit.value()) > int(self.max_unit.value()) else int(self.min_unit.value())
    		#self.min_unit.setValue(1) if (int(self.min_unit.value()) < 6) and (int(self.min_unit.value()) !=int(self.max_unit.value()))  else int(self.min_unit.value())
    		self.min_unit.setValue(int(self.max_unit.value())) if int(self.min_unit.value()) >  int(self.max_unit.value())\
    		                                                   else int(self.min_unit.value())
    	# check if DNA or RNA sequence for Default (fast, less accurate) engine selection
    	if ((self.DNA_button.isChecked() == True  or self.RNA_button.isChecked() == True ) and self.default_eng.isChecked() == True):
    		self.max_unit.setValue(10) if int(self.max_unit.value()) > 10 else int(self.max_unit.value())
    		self.min_unit.setValue(self.max_unit.value()) if int(self.min_unit.value()) > int(self.max_unit.value()) else int(self.min_unit.value())
    		#self.min_unit.setValue(1) if (int(self.min_unit.value()) < 10) and (int(self.min_unit.value()) !=int(self.max_unit.value()))  else int(self.min_unit.value())
    		self.min_unit.setValue(int(self.max_unit.value())) if int(self.min_unit.value()) >  int(self.max_unit.value()) \
    		                                                      else int(self.min_unit.value())
    	# check if protein sequence for default( less accurate but fast search engine selection
    	if ((self.Protein_button.isChecked() == True ) and self.default_eng.isChecked() == True):
    		self.max_unit.setValue(3) if int(self.max_unit.value()) > 3 else int(self.max_unit.value())
    		self.min_unit.setValue(self.max_unit.value()) if int(self.min_unit.value()) > 3 else int(self.min_unit.value())
    		#self.min_unit.setValue(1) if (int(self.min_unit.value()) < 3) and (int(self.min_unit.value()) !=int(self.max_unit.value()))  else int(self.min_unit.value())
    		self.min_unit.setValue(int(self.max_unit.value())) if int(self.min_unit.value()) >  int(self.max_unit.value()) \
    		                                                      else int(self.min_unit.value())
    	# check if protein sequence for Accurate engine selection
    	if ((self.Protein_button.isChecked() == True ) and self.accurate_eng.isChecked() == True):
    		self.max_unit.setValue(4) if int(self.max_unit.value()) > 4 else int(self.max_unit.value())
    		self.min_unit.setValue(self.max_unit.value()) if int(self.min_unit.value()) > 4 else int(self.min_unit.value())
    		self.min_unit.setValue(1) if (int(self.min_unit.value()) < 4) and (int(self.min_unit.value()) !=int(self.max_unit.value()))  else int(self.min_unit.value())
    		self.min_unit.setValue(int(self.max_unit.value())) if int(self.min_unit.value()) >  int(self.max_unit.value()) \
    		                                                      else int(self.min_unit.value())

    	# check if DNA or RNA sequence for Default (fast, less accurate) engine selection
    	if self.dict_mode.isChecked() == True:
    		self.min_unit.setValue(self.max_unit.value()) if (int(self.min_unit.value()) < 30) and (int(self.min_unit.value()) !=int(self.max_unit.value())) \
                                                                  else int(self.min_unit.value())
    		self.max_unit.setValue(30) if int(self.max_unit.value()) > 30 else int(self.max_unit.value())
    		self.min_unit.setValue(self.max_unit.value()) if int(self.min_unit.value()) > 30 else int(self.min_unit.value())
    		self.min_unit.setValue(self.max_unit.value()) if (int(self.min_unit.value()) > int(self.max_unit.value())) and (int(self.min_unit.value()) !=int(self.max_unit.value()))  else int(self.min_unit.value())
    		self.min_unit.setValue(int(self.max_unit.value())) if int(self.min_unit.value()) >  int(self.max_unit.value()) \
                                                                  else int(self.min_unit.value())

    #---------------------------all checks in the DNA/RNA/Protein ticks done
    def Activate_flanking_seq(self):                                           

        if (self.DNA_button.isChecked() == False) and (self.PrimerFlag.isChecked() == True):
            QMessageBox.about(self, "OSTRFPD Warning!", "The make primer option is only valid for DNA sequence source")
            self.PrimerFlag.setChecked(False)

        if (self.PrimerFlag.isChecked() == True):
            self.left_flanking_flag.setChecked(True)
            self.right_flanking_flag.setChecked(True)
            self.left_flank_seq.setValue(85) if (int(self.left_flank_seq.value())   <=0) else None
            self.right_flank_seq.setValue(85) if (int(self.right_flank_seq.value()) <=0) else None

    def open_folder(self):
        curfolder =  os.curdir   #  get teh current directory 
        warnings.filterwarnings("ignore") 
        os.system('xdg-open "%s"' % curfolder) if os.name != "nt" else os.startfile(curfolder)  # opeen in explorer                                                                         #current folder for linux or windows
        

    def GUI_module(self):                                  # run the microsatelliet from GUI module
        global fsc_string                                  # for the fsc_string global value change
        global primer_saveFileName                         # For the global primer save file pathname
        global dict_file                                   # makes a list from dictionary based on length         
        
        imperfect_only = False

        self.log_window.append("\n" + str(datetime.now()) + "\n[Scanning task started]")  # start task log window
        os.system('cls' if os.name=='nt' else 'clear') if \
                 (self.console_screen_flag.isChecked()== True) else  None    # system independent console screen clear ( only works in run mode not debug mode)

        if (self.PrimerFlag.isChecked() == True):
            self.left_flanking_flag.setChecked(True)
            self.right_flanking_flag.setChecked(True)

                                                            # remove the file:/// term for windows os 
        targetfile = str(self.input_file.toPlainText()) if os.name != "nt" else \
                                                           str(self.input_file.toPlainText()).strip().replace("file:///", "")
        outfile =str(self.output_file.toPlainText()) if os.name != "nt" else \
                                                          str(self.output_file.toPlainText()).strip().replace("file:///", "")

        if (self.save_primer.toPlainText().strip() == "" and self.PrimerFlag.isChecked() == True):
            QMessageBox.about(self, "OSTRFPD Alert!", "Please enter a valid primer_result file name [e.g: primer_results.txt]")            # GUI command for finish task
            self.save_primer.setText ("primer_result.txt")
            return

        if (self.min_fix.value() == 0 and self.dict_mode.isChecked() == True):
            QMessageBox.about(self, "OSTRFPD Alert!", "Missing Input! Please set the value for Minimum Repeat Number to use with Dictionary search.")
            return

        self.misa_string.setText("") if ( (self.dict_mode.isChecked() == True) or ( self.min_fix.value() > 0 ) ) else None
        self.check_motif()
        primer_saveFileName = self.save_primer.toPlainText().strip()
        unit_motif_len = int(self.max_unit.value())
        unit_motif_min_len = int(self.min_unit.value())
        min_repeat= int(self.min_fix.value())
        uFix_flag = True if int(self.min_unit.value()) == int(self.max_unit.value()) else False
        imperfect_repeate = int(self.max_imperfect.value())
        left_ext_len = int(self.left_flank_seq.value())   if ( self.left_flanking_flag.isChecked() == True ) else 0
        right_ext_len = int(self.right_flank_seq.value()) if ( self.right_flanking_flag.isChecked() == True ) else 0
        im_occur_penalty =int(self.missmatch_open.value())
        miss_occur_penalty = int(self.non_seed_strat_penalty.value())
        misa_string =  str(self.misa_string.toPlainText().strip())
        mmp_val= int(self.missmatch_extension_penalty.value())
        strictFlag = True if  (self.strict_flag.isChecked() == True)   else False
        if (self.DNA_button.isChecked() == True):
            dnaFlag = "dna"
        if (self.RNA_button.isChecked() == True):
            dnaFlag = "rna"
        if (self.Protein_button.isChecked() == True):
            dnaFlag = "protein"
        engFlag        = True if (self.accurate_eng.isChecked()== True) else False

        if (self.stat_flag.isChecked()  == False): # opposit login cuz true is ON by default
        	stats_flag     = True
        	report_flag    = True
        else:
        	stats_flag     = False
        	report_flag    = False    

        if (unit_motif_min_len  > unit_motif_len):                                          # Set min unit motif length = unit motif length if greater than maximum
             unit_motif_min_len = unit_motif_len

        show_scoreflag = True if (self.sshow_flag.isChecked()  == True) else False
        fsc_string     = str(self.fsc_string.toPlainText().strip())
        fsc_flag       = True if str(self.fsc_string.toPlainText().strip()) != "" else False
        sdout_flag     = True if (self.show_console.isChecked() == False) else False        # opposite login cuz true is ON by default
        min_gap        = int(self.min_gap.value())
        gapflag        == True if int(self.min_gap.value()) != 0 else False
        revComp        = True if self.rev_comp_flag.isChecked() == True else False    
        align_flag     = True if (self.output_format.currentItem().text()) == "Alignment" else False
        if (self.output_format.currentItem().text()) == "Imperfect Alignment":
            align_flag     = True
            imperfect_only = True
        fasta_flag     = True if (self.output_format.currentItem().text()) == "FASTA" else False
        sim_value      = int(self.min_perfection.value()) 
        if self.dict_mode.isChecked() == True:
            dict_file      = str(self.dict_string.toPlainText()) if os.name !='nt' else str(self.dict_string.toPlainText()).strip().replace("file:///", "")
        else:
            dict_file      = ""   
        stdFlag        = True if (self.motif_type.currentItem().text()) == "Full:[rev.complement + cyclic]" else False 
        write_header   = True  # for first
        make_primer    = True if (self.PrimerFlag.isChecked()== True) else False
        l_tag          = self.lTag_seq.currentItem().text() if (self.lTag_flag.isChecked() == True) else "" 
        r_tag          = self.rTag_seq.currentItem().text() if (self.rTag_flag.isChecked() == True) else "" 

        # Primer parameters setting ------------------------------------------------------------------------------
        prm_range    = str(self.prm_range.toPlainText())
        prm_opt_size = str(self.prm_opt_size.toPlainText()) 
        prm_min_size = str(self.prm_min_size.toPlainText()) 
        prm_max_size = str(self.prm_max_size.toPlainText()) 
        prm_opt_tm   = str(self.prm_opt_tm.toPlainText()) 
        prm_min_tm   = str(self.prm_min_tm.toPlainText()) 
        prm_max_tm   = str(self.prm_max_tm.toPlainText()) 
        prm_tm_diff  = str(self.prm_tm_diff.toPlainText()) 
        prm_opt_gc   = str(self.prm_opt_gc.toPlainText()) 
        prm_min_gc   = str(self.prm_min_gc.toPlainText()) 
        prm_max_gc   = str(self.prm_max_gc.toPlainText())
        prm_poly     = str(self.prm_poly.toPlainText()) 
        #----------------------------------------------------------------------------------------------------------
        if (imperfect_only == True and imperfect_repeate == 0):
        	QMessageBox.about(self, "OSTRFPD Alert!", "Please set a Non-zero positive value for the imperfect repeat. Normally length of motilf * 2 ")            # GUI command for finish task
        	return

       #----------------------------------------------------------------------------------------------------------
        retval =Search_Microsatellite (targetfile, outfile, unit_motif_len, unit_motif_min_len, min_repeat, uFix_flag, imperfect_repeate
                                              , left_ext_len, right_ext_len, im_occur_penalty, miss_occur_penalty 
                                              , misa_string, mmp_val, strictFlag, dnaFlag, engFlag, stats_flag 
                                              , show_scoreflag, fsc_flag, fsc_string, fasta_flag, sdout_flag, gapflag
                                              , min_gap, revComp, report_flag, align_flag, sim_value
                                              , dict_file, stdFlag, write_header, make_primer, l_tag, r_tag, imperfect_only
                                              , prm_range, prm_opt_size, prm_min_size, prm_max_size, prm_opt_tm, prm_min_tm, prm_max_tm
                                              , prm_tm_diff, prm_opt_gc, prm_min_gc, prm_max_gc, prm_poly)    
      
        self.log_window.append("\n" + str(datetime.now()) + "\n[Scanning task finished!]")
        QMessageBox.about(self, "OSTRFPD Alert!", "Sequence scanning completed!")            # GUI command for finish task


        #UNCOMMENT To Draw motif frequency chart if search results are not None
        #motif_chart(retval) if retval != [] else None

        return 

#-------------------Function gets the alignment using BioPython global alignment for the segment of microsatellite 
def get_alignemnt(seqID, mss, mrn, ms_type, msp, mep, base_score, non_seed_missmatch, seed_missmatch, LFseq, RFseq ):
    diff = int((len(mss) - mss.count(ms_type)* len(ms_type))/len(ms_type)) # predict minimum abount of extra motif to be added in perfect sequence as equivalent
    refSeq = (mss.count (ms_type) + diff ) * ms_type                       # same strategy
    alignments = pairwise2.align.localmx(mss, refSeq, 2, -1)

    al1,al2, score, begin, end = alignments [0][1],alignments [0][1],alignments [0][2],alignments [0][3],alignments [0][4]
    consensus, match , result =[] , [], ""
    # general dictionaly to state the strength of asociation
    cons_dic = {'CG': 'S', 'AT': 'W', 'AC': 'M', 'AG': 'R'}  if dnaFlag == "dna" else None         
    cons_dic = {'CG': 'S', 'AU': 'W', 'AC': 'M', 'AG': 'R'}  if dnaFlag == "rna" else None
    cons_dic = {'KG': 'S', 'GI': 'W', 'LL': 'M', 'AH': 'R'}  if dnaFlag == "protein" else None

    for a, b in zip(alignments[0][0],alignments[0][1]):
        if a == b:
            consensus.append(a)
            match.append('|')
        elif a == '-':
            consensus.append(b)
            match.append(' ')
        elif b == '-':
            consensus.append(a)
            match.append(' ')
        else:
            key = ''.join(sorted(a + b))
            consensus.append(cons_dic[key])

    m="".join(match)
    s=[]
    if len(LFseq) >10 and len(RFseq)>10:                              # show 10 bp flanking sequence of possible
        s.append ("[Identified repeat]" +'\n')
        s.append( LFseq[len(LFseq)-10::]+ "==="+ alignments[0][0]+ "===" + RFseq[:10] + '\n' )
        s.append(len(LFseq[len(LFseq)-10::]+ "===") * " "  + m+'\n')                        # manage space 
        s.append( LFseq[len(LFseq)-10::]+ "==="+ alignments[0][1] +"===" + RFseq[:10] + '\n')
        s.append("[~Equivalent perfect repeat]")
    else:
        s.append ("[Identified, computed repeat]" +'\n')
        s.append( LFseq[len(LFseq)-10::]+ "==="+ alignments[0][0]+ "===" + RFseq[:10] + '\n')
        s.append(len(LFseq[len(LFseq)-10::]+ "===") * " "  + m+'\n')                        # manage space 
        s.append( LFseq[len(LFseq)-10::]+ "==="+ alignments[0][1] +"===" + RFseq[:10]     )
        s.append("\n[~Equivalent perfect repeat]")      

    alignedSeqs="".join(s)
    result =result + "\n" +  ('Sequence ID            :' + str(seqID))
    result =result + "\n" +  ('Motif type             :' + str(ms_type))
    result =result + "\n" +  ('Start position (bp|aa) :' + str(msp) )
    result =result + "\n" +  ('End position   (bp|aa) :' + str(mep))
    result =result + "\n" +  ('Repeat number          :' + str(mrn))
    result =result + "\n" +  ('Repeat length          :' + str(mep-msp))
    if base_score  >=0:                                               # should write like this if not all above will be neglected
        result = result + "\n" + ('Repeat type            :Perfect') 
    elif base_score <0:
        result = result + "\n" + ('Repeat type            :Imperfect')
    result =result + "\n" +  ('Local alignment score   :' + str(score))
    result =result + "\n" +  ('OSTRFPD custom score    :' + str(base_score))
    result =result + "\n" +  ('Non-motif missmatch     :' + str(non_seed_missmatch))
    result =result + "\n" +  ('Motif missmatch         :' + str(seed_missmatch))
    result =result + "\n" +  ('Computed consensus      :' + ''.join(consensus))
    result =result + "\n" +  ('\n')
    result =result + "\n" +  (alignedSeqs)
    result =result + "\n" +  ("\n")
    result =result + "\n" +  ('Observed  : ' + alignments[0][0]  )
    result =result + "\n" +  ('Perfect   : ' + alignments[0][1]  )
    result =result + "\n" +  ('Consensus : ' + ''.join(consensus))
    result =result + "\n" + "__________________________________________________________"
    return (result +"\n")

#------------------- Primer3_configuration file for primer3_core
def primer3_config(   primer_sequence_id,sequence,target,primer_product_size_range, primer_opt_size, primer_min_size
                    , primer_max_size, primer_opt_tm, primer_min_tm, primer_max_tm, primer_max_diff_tm
                    , primer_opt_gc_percent, primer_min_gc_percent, primer_max_gc_percent, primer_max_poly_x
                    , primer_num_return):

    primer3_string ="PRIMER_SEQUENCE_ID=" + primer_sequence_id 
    primer3_string +="\n" + "SEQUENCE="   + sequence
    primer3_string +="\n" + "TARGET="     + target 
    primer3_string +="\n" + "PRIMER_PRODUCT_SIZE_RANGE=" + primer_product_size_range 
    primer3_string +="\n" + "PRIMER_OPT_SIZE=" + primer_opt_size
    primer3_string +="\n" + "PRIMER_MIN_SIZE=" + primer_min_size
    primer3_string +="\n" + "PRIMER_MAX_SIZE=" + primer_max_size
    primer3_string +="\n" + "PRIMER_OPT_TM="   + primer_opt_tm
    primer3_string +="\n" + "PRIMER_MIN_TM="   + primer_min_tm
    primer3_string +="\n" + "PRIMER_MAX_TM="         +    primer_max_tm
    primer3_string +="\n" + "PRIMER_MAX_DIFF_TM="    +    primer_max_diff_tm
    primer3_string +="\n" + "PRIMER_OPT_GC_PERCENT=" +    primer_opt_gc_percent
    primer3_string +="\n" + "PRIMER_MIN_GC_PERCENT=" +    primer_min_gc_percent
    primer3_string +="\n" + "PRIMER_MAX_GC_PERCENT=" +    primer_max_gc_percent
    primer3_string +="\n" + "PRIMER_MAX_POLY_X="     +    primer_max_poly_x
    primer3_string +="\n" + "PRIMER_NUM_RETURN="     +    primer_num_return
    primer3_string +="\n" + "PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION="     +    target
    primer3_string +="\n" + "="

    return primer3_string

#------------------------------------------function to create or add the next primer using primer3_core.exe (same foelder as python file)
def write_primer_seq(primer3_file, config_data, sequence, write_header, l_tag, r_tag, ms_type, mrn):

    if (write_header == True):                                              # write the header on new priemr list file
        with open(primer3_file,"w") as f:
            f.write ("# OMTRFPD v0.01 Automated primer design output. DateTime stamp: " + str(datetime.now()) + "\n" )
            f.write ("sequence_id" + "\t" + "left_primer" + "\t" +  "right_primer" 
                                   + "\t" + "left_Tm"     + "\t" +  "right_Tm" + "\t" + "size(bp)"
                                   + "\t" + "microsatellite_motif" + "\t" + "repeat_number" + "\t" + "full_sequence" )    

    with open("tmp_config","w") as f:                                                        # write the configuration file 
        f.write(config_data) 

    batch_command= ("primer3_core.exe < tmp_config"  ) if os.name == "nt" else \
                   ("./primer3_core   < tmp_config"  )                                       # batch command for windows and linux for os dependent binaries

    try:
        result = str(subprocess.check_output(batch_command, shell=True))                     # convert to string for the byte array type return

    except:
        return "-1" 																	     # is something wrong and cannot make primer return -1
	# solution to OS-dependent and version related issues with file formatting and primer3_core                                                        
    if os.name == "nt":
    	ans   = result.split("\r\n") if (version_info[0] < 3) else (result.split("\\r\\n"))  # for python 2 and 3 both compability due to return type from sub-process in linux                                                               # for linux version    						                                                                    
    else:
    	ans   = result.split ("\\n")                                                         # for linux version  or otehrs (not tested)              	
    #=====================
    if result.count("PRIMER_PRODUCT_SIZE=") == 0:                                            # so no primer is being formed..so exit 
        return "-1"
    lst =[]
    for x in ans:
        if x.count("PRIMER_SEQUENCE_ID") ==1:           
            lst.append((x.strip()).split("=")[1])                                       # split the value ofrm '=' and get the right side value
        if x.count("PRIMER_LEFT_")  ==1  and (x.count("_SEQUENCE=") == 1):
            lst.append(l_tag + (x.strip()).split("=")[1])                               # get priemr sequence and add left primer tags if given
        if x.count("PRIMER_RIGHT_") ==1 and (x.count("_SEQUENCE=") == 1):
            lst.append(r_tag + (x.strip()).split("=")[1])                               # get priemr sequence and add right primer tags if given
        if x.count("PRIMER_LEFT_") ==  1 and x.count("_TM=") ==1:
            lst.append((x.strip()).split("=")[1]) 
        if x.count("PRIMER_RIGHT_") and x.count("_TM=") == 1:
            lst.append((x.strip()).split("=")[1])               
        if x.count("PRIMER_PRODUCT_SIZE") and x.count ("PRIMER_PRODUCT_SIZE_RANGE") ==0:
            lst.append((x.strip()).split("=")[1])                                       # put the new line carriage after each end set
    

    # write each primer outputinformation 
    with open(primer3_file,"a+") as f:
        lst.append (str(ms_type))                                                       # append standarized microsatelliet motif
        lst.append (str(mrn))                                                           # append microsatellite repeat number
        lst.append (sequence)                                                           # append full sequence  (LF+ mms + RF)
        f.write ("\n" + "\t".join(lst))                                                 # join and write the primer sequence.

    return lst

#-------------------------------Makes the motif sequence list from the custom sequence supplied in dictionary (fixed length)
def dict_pat(uml):
    global dict_file                                                                    

    with open(dict_file) as f:
        content = f.readlines()
        content = [x.strip() for x in content] 
    lst=[]
    for x in content:
        lst.append(x)  if len(x) == uml  else None
    return lst

#-------------------------------Getss the report 
def get_report(text_string, dnaFlag):
    frequency = {}
    match_pattern = re.findall(r'\b[A-Z]{1,}\b', text_string)
    for word in match_pattern:
        count = frequency.get(word,0)
        frequency[word] = count + 1
    frequency_list = frequency.keys() 
    report_list=[]
    for words in frequency_list:
    	mem_word =words
    	words = (words +" [" + seq3(words) +"]") if (dnaFlag == "protein") else words
    	report_list.append( [words, frequency[mem_word]] )
    report_list= sorted(report_list, reverse =True)
    return report_list

#-----------------------------Detect string rotation simplest way to drop down the pattern list
def is_rotation(str1, str2):  
	if (len(str1) == len(str2) and str1 in str2*2):
		return True
	else:
		return False 

#------------------------------Detect if the word is multiple of itself and returns True is yes
def repeats(string):                                          
    for x in range(1, len(string)):
        substring = string[:x]
        if substring * (len(string)//len(substring))+(substring[:len(string)%len(substring)]) == string:
            return True

#-------------------------------Check if the motif is palindrome or Not
def isPalindrome(motif_type):                                  
    for i in range(0, int(round(len(motif_type)/2))):          #  have to use int(round(float())) go get the integer becuase we impotred devision
        if motif_type[i] != motif_type[len(motif_type)-i-1]:
            return False
    return True
 
#-------------------------------#Returns the best alphabetically arranged Motif seed of all possible cyclic rotation 
def get_std_seed(seed):                                        
    a = seed
    b = len(a)
    tmp =[]
    for i in range (b):
        c = a[i:]+a[:i]
        tmp.append(c)
        tmp = sorted(tmp)
    return tmp[0] if tmp else None

#------------------------------#Standarize the motifs so that both original and reverse compelment will be noted as same in output
def std_motif(ms_type):      
    lst = []
    lst.append(ms_type)
    dna = Seq(ms_type)
    lst.append(str(dna.reverse_complement()))
    lst =sorted(lst)
    return lst[0]

#-------------------------------Combine the tr_find function with flag to compute if the flanking sequence has repeat or not                    
def flanking_check(flankSeq):                                 
    global tr_flag
    tr_flag =False                                            # reset/ or set the global flag to False (default) before checking
    check = tr_find(flankSeq)                                 # cheanges the fsc_flag value to true if found
    return tr_flag                                            # recursive function to check all the repeats and return boolean

#-------------------------------Recursive function for (suitable for short few kb sequences) to find tandem repeats of minimum unit length (fsc_uml) and minimum repeat (fsc_umr)
def tr_find(seq,past_pos=0, result=[]):                       
   global fsc_uml
   global fsc_umr  
   global tr_flag                                             # must be shared between both main module and  each sub module
   global fsc_string
   
   if (fsc_string.strip() == ""):
        return False

   tmp = fsc_string.strip().split(",")                        # in smeewhat MISA format x;y
   fsc_uml = int(tmp[0])                                      # gets minimum motif length
   fsc_umr = int(tmp[1])                                      # gets minimum motif reepat 

   tr_found=re.search(r'(.+)\1+',seq)                         # search for all form of repats in flanking sequence 
   if tr_found:
      tr_strt,tr_ends=tr_found.span()
      sub=seq[tr_strt:tr_ends]
      ind = (sub+sub).find(sub, 1)
      sub=sub[:ind]
      rn=((tr_ends-tr_strt)/len(sub))
      #============================
      if (len(sub) >= fsc_uml) and (seq.count(fsc_umr * sub) >= fsc_umr):  # for avoiding any subtext match et ATTTTATTTTATTTT for fsc =1;4 but regx fails to detect
         tr_flag = True
         return tr_find(seq[tr_ends:],past_pos)
      #============================                            # extra layer of protection for loger reepats that may not be detect by count fut will be by regx
      if (len(sub)>= fsc_uml) and (rn >= fsc_umr):             # minimum motif length and repeat condition to be True
        result.append([sub,rn,(tr_strt + past_pos + 1,tr_ends + past_pos), seq[tr_strt:tr_ends]])
        tr_flag = True
      past_pos+=tr_ends
      return tr_find(seq[tr_ends:],past_pos)
   else:
      return result

#-----------------------------------Get the minimum motif in MISA format and returns an equivalent array
def getMISA(misa_string):                                    
    misa_string =misa_string.split(",")                      # split based on the delimiter ',' and make array
    count =0
    misa_motif=[]
    for x in misa_string:
        misa_motif.append(int(x))

    return misa_motif                                        

#-----------------------------------Generates reverse complements using biopython function
def seq_revComp(sequence):                                   
    dna = Seq(sequence)
    return str(dna.reverse_complement())

#----------------------------------Generates complements sequence using biopython
def seq_Comp(sequence):                                      
    dna = Seq(sequence)
    return str(dna.complement())
                                                             
#----------------------------------Generate dictionary for [uml sized]  length motif  without any repating or rotating sequences (unique) for given length
def patternList(me_block, unit_motif, engFlag):              
    global dict_file                                         # makes a list from dictionary based on length 
    char_set = me_block
    complete_list = []
                                                             
    if (engFlag == True) and (dict_file == ""):              # DICTIONAY GENERATION. use slow but highly efficent built-in pattern generator (brute force)
        for current in range(unit_motif-1,unit_motif):       # Gets the all possible set for each unit motif length 
            a = [i for i in char_set]
            for y in range(current):
                a = [x+i for i in char_set for x in a]
            complete_list = complete_list+a

    elif (engFlag == False) and (dict_file !=""):
        complete_list = dict_pat(unit_motif)                 # read the motif list from the dictionay file if dictionay file option is supplied      

    else:                                                    # fast but low efficent all character set..less relable for longer patterns   
        tmp_motif_set=cwr(me_block, unit_motif )
        for x in tmp_motif_set:
            motif= "".join(string for string in x)
            complete_list.append(motif)

    count,flist =0, []
    for x in range(len(complete_list)):
        myflag =0
        for y in range(x+1, len(complete_list)-1):
            if (is_rotation(complete_list[x], complete_list[y]) == True):               
                myflag =1
            if (len(complete_list[x]) >1) and (isPalindrome(complete_list[x]) == True):
                myflag =1

        if (myflag ==0) and (repeats(complete_list[x]) != True):        # Flag should be zero and remain should not be either palendromic or repetative avoid tautoligical rotatory term motifs AAT TAA and palendromic motifs as ATA AATAA
            if len(set(complete_list[x])) >1 and (len(complete_list[x]) >1):
                flist.append(complete_list[x])
            elif len(set(complete_list[x])) >1 and (len(complete_list[x]) == 1):
                flist.append(complete_list[x])
            elif len(complete_list[x]) == 1:
                flist.append(complete_list[x])

    return flist                                                                     # retursn as list one dimension

#-----------------------------Builds the regular expression pattern searech of uml length motif only
def build_regex(uml,min_repetitions,imperfect, me_block, engFlag):                                 
    plist = patternList(me_block, uml, engFlag )                                     # run-time generate the pattern list for the search
    pattern = "|".join("(?:" + "".join(i) + "){" + str(min_repetitions if (uml == 1 and imperfect <= 0) 
        else min_repetitions - 1) + ",}" for i in plist)
    #print("\nTest pattern: ",  pattern)                                              # Uncomment to print out each search pattern                         
    return pattern


#-----------------------------Includes the list of Imperfect repeats
def include_imperfect_repeats(rslts, imperfect, strict, seq, gapflag, min_gap):
    if len(rslts) == 0:
        return rslts

    rslts = sorted(rslts, key=lambda x: x[1])
    _rslts = [rslts[0]]

    for i in range(1, len(rslts)):

        # Check for teh seed and the distance between the two repeats
        if (_rslts[-1][3] == rslts[i][3]) and (rslts[i][1] - _rslts[-1][2] <= imperfect):
            if strict==True:
                # Get the sequence between the two repeats
                sub_sequence = seq[_rslts[-1][2]: rslts[i][1]]
                if len(set(sub_sequence) - set(rslts[i][3])) > 0:          # for checking only the nucleotide if present in seed
                    continue
                # Extend the end
                _rslts[-1][2] = rslts[i][2]
                # Extend the sequence
                if imperfect == 0:
                    _rslts[-1][4] += rslts[i][4]
                
                elif strict==False:
                    sub_sequence = seq[_rslts[-1][2]: rslts[i][1]]
                    _rslts[-1][2] = rslts[i][2]

                else:
                     _rslts[-1][4] = seq[_rslts[-1][1]: _rslts[-1][2]]
            
            if strict==False:
                # Get the sequence between the two repeats
                sub_sequence = seq[_rslts[-1][2]: rslts[i][1]]
                _rslts[-1][2] = rslts[i][2]
                # Extend the sequence
                if imperfect == 0:
                    _rslts[-1][4] += rslts[i][4]
           
                else:
                     _rslts[-1][4] = seq[_rslts[-1][1]: _rslts[-1][2]]    
        else:
            # Append a new one
            if gapflag == True:               # add only if minimum gap exists flag is true and satisfies min gap
                _rslts.append(rslts[i]) if (rslts[i][1] - _rslts[-1][2]) >= min_gap else None
            else:
                _rslts.append(rslts[i])   # just not add if not satisfy gap condition ( but may satisfy imperfect condition)
    return _rslts

#-----------------------------Returns the penalty scores for imperfect repeats
def missmatch_penalty(mss, ms_type, me_block, mmp_val, im_occur_penalty, miss_occur_penalty):    # returns the gap penalty of non_seed item mismatcg

    non_seed = set(me_block) -set(ms_type)                                                    # entire unique nucleotide/ or protein set - set items in seed
    non_seed_count=0

    for x in non_seed:
        non_seed_count=non_seed_count + mss.count(x)

    non_seed_score = non_seed_count * mmp_val                                            # get the non_seed missmatch score only 

    if mss in ( (mss.count(ms_type)+ 2)  * ms_type ):                                   # strategy to check is mss is perfectly inside the length pf mss + 2* unit_motif   to get edge correction
        base_score =0
    else:
        base_score = len(mss) - (mss.count(ms_type) * len(ms_type) )                    # get basic score [ microsatellit length - length of perfect match within microsatellite] to detect imperfectness
    
    seq_identity = round( (100*(len(mss)-base_score) /len(mss)), 3 )                    # rough percentage of sequence elements (nucleotide or proteins) involved in perfect repetative sequence for total length of given mss 
    seed_missmatch_count =  base_score - non_seed_count                                 # get total missmatches type contaning witin the seed (rounded to 2 decimal palce)

    if   (base_score != 0) and (im_occur_penalty != 0 ) and (miss_occur_penalty == 0):
        score = base_score + im_occur_penalty
    elif (base_score != 0) and (im_occur_penalty != 0) and  (miss_occur_penalty != 0):
        score = base_score + im_occur_penalty + miss_occur_penalty + non_seed_score     # by definition im_occur_penalty wwill be added by default for non_seed missmatch penalty and opening of non_seed missmatcg
    else:
        score= base_score

    return [score, base_score, non_seed_count, seed_missmatch_count, seq_identity]                                          # returns total score, base score, non_seed_count as a list

#----------------------Main search modulle-----------------------------------------------------------
def Search_Microsatellite(targetFile, outfile, unit_motif_len, unit_motif_min_len, min_repeat, uFix_flag, imperfect_repeate
                                    , left_ext_len, right_ext_len, im_occur_penalty, miss_occur_penalty
                                    , misa_string, mmp_val, strictFlag, dnaFlag, engFlag, stats_flag 
                                    , show_scoreflag, fsc_flag, fsc_string, fasta_flag, sdout_flag
                                    , gapflag, min_gap, revComp, report_flag, align_flag, sim_value
                                    , dict_file, stdFlag, write_header, make_primer, l_tag, r_tag, imperfect_only
                                    , prm_range, prm_opt_size, prm_min_size, prm_max_size, prm_opt_tm, prm_min_tm
                                    , prm_max_tm,prm_tm_diff, prm_opt_gc, prm_min_gc, prm_max_gc, prm_poly):

    
    #--------------------Major variables declerations
    global primer_saveFileName      # for the global primer save file name
    global dict_warning

    start_time = time.time()
    if dnaFlag == "dna":                                                         # for the DNA sequence
        me_block = "ATGC"
    elif dnaFlag == "rna":                                                       # for the RNA sequence
        me_block = "AUGC" 
    else:
        me_block = "ACDEFGHIKLMNPQRSTUVWY"                                       # for the protein sequence

    uml                     = unit_motif_len
    min_repetitions         = min_repeat
    imperfect               = imperfect_repeate
    motif_log               = ""                                                # for logging motif
    total_ms_count          = 0                                                 # total ms count for log
    total_ms_coverage       = 0                                                 # get total coverage
    total_seq_scan_length   = 0                                                 # total sequence scanned length
    total_fasta_seq         = 0                                                 # get total number of fasta file (in multi fasta)
    motif_freq_byLength     = array.array( "i", (0 for i in range (0,uml+2)))   # set the array of intiger from 1 to  uml+2 (for safety) initalized as zero 
    motif_coverage_byLength = array.array( "i", (0 for i in range (0,uml+2)))   # get coverage length of each motif
    primer_count            = 0 
    org_motif               = ""                                                # maintain org (rotatory config) of last motif (for alignemnt mode)
    report_data             = "" 
    GC_avg                  = 0.00
    imperfect_count         = 0                                                 # count the imperfect repeat
    ms_gc                   = 0

    #-----------------------------------------------------
    if uFix_flag == True:                                                       # set up the fix length of motif to extract
        uFix = uml
    else:
        uFix = unit_motif_min_len                                               # for unit motif to uml (min to uml)

    with open(outfile,"w") as f:

        if os.path.splitext(targetFile)[1] == '.gz':              # Handel gunzipped file
            handel =gzip.open(targetFile, "rt") if (version_info[0] > 2   ) else \
                                                gzip.open(targetFile, "rb")      # Solved Fix gzip issuesh (python 3.2 and 2.x) 
        else:
            os.path.splitext(targetFile)[1] == '.fasta'
            handel = targetFile
        #------------------------------------------------------------------------------# Start of pre-scan messages 
        # some header information comment on file and input and target file info 
        built_info = (25* "__", "\nStarting new task\n")
        built_info = "#OSTRFPD v0.01 generated output\n#Scan output [datetime]  : " + str(datetime.now())   
        built_info = (built_info + "\n") +( "#Processing file         : [%s]"% targetFile)
        built_info = (built_info + "\n") +( "#Output result file      : [%s]"% outfile)

        misa =[]                                                                       # Initalize MISA array
        if (misa_string == ""):                                                        # default MISA vlaues for DNA nucleotide
            if dnaFlag == "dna" or (dnaFlag =="rna"):
                misa_string = "14,7,5,4" + uml * ",4"                                  # put some extra long MISA patter just to be secure of length limit      
            else:
                misa_string = "5,4,3,3"  + uml * ",3"                                  # put some extra motif length just to make sure it covers 1...uml
        
        if uFix_flag == True:
            built_info = (built_info + "\n") + ("#Motif search range fixed: [%d]"%uml)
        else:
            built_info = (built_info + "\n") + ("#Motif search range      : ["+str(uFix)+ "-" + str(uml)+ "]")

        if dnaFlag   =="dna":
            built_info = (built_info + "\n") + ("#Processing sequence type: [DNA]")
        elif dnaFlag =="rna":
            built_info = (built_info + "\n") + ("#Processing sequence type: [RNA]")
        else:
            built_info = (built_info + "\n") + ("#Processing sequence type: [Protein]")

        if imperfect_repeate <= 0:
            built_info = (built_info + "\n") + ("#Search pattern Mode     : [Perfect]")
        else:
            built_info = (built_info + "\n") + ("#Search pattern Mode     : [Imperfect]")

        if fasta_flag == True:
            built_info = (built_info + "\n") + ("#Output file format type : [Fasta]")
        elif align_flag == True and imperfect_only == False:
            built_info = (built_info + "\n") + ("#Output file format type : [Alignment]")
        elif align_flag == True and imperfect_only == True:
            built_info = (built_info + "\n") + ("#Output file format type : [Imperfect alignment only]")
        else:
            built_info = (built_info + "\n") + ("#Output file format type : [Custom]")

        if (report_flag == True) and (fasta_flag == False) :
            built_info = (built_info + "\n") + ("#Attach report to output : [Yes]")
        elif (report_flag == True) and (fasta_flag == True):
            built_info = (built_info + "\n") + ("#Attach report to output : [Not applicable]")
        elif (report_flag == False):
            built_info = (built_info + "\n") + ("#Attach header and report: [No]")

        if revComp ==True and ( (fasta_flag == False and align_flag == False ) and  (dnaFlag == "dna" or dnaFlag == "rna") ) :
            built_info = (built_info + "\n") + ("#Reverse strand sequence : [Included]")
        elif revComp ==True and (fasta_flag == True or align_flag == True):
            built_info = (built_info + "\n") + ("#Reverse strand sequence : [Not applicable]")
        else:
            built_info = (built_info + "\n") + ("#Reverse strand sequence : [Not included]")

        built_info = (built_info + "\n")     + ("#Similarity (Percentage) : [%d Minimum]"%sim_value)

        if stats_flag ==True:
            built_info = (built_info + "\n") + ("#Show basic statistics   : [Yes]")
        else:
            built_info = (built_info + "\n") + ("#Show basic statistics   : [None]")

        if sdout_flag ==True:
            built_info = (built_info + "\n") + ("#Standard console output : [Activated]")
        else:
            built_info = (built_info + "\n") + ("#Standard console output : [Forced OFF]")

        if fsc_flag ==True:
            built_info = (built_info + "\n") + ("#Flanking sequence filter: [Activated: " + fsc_string + "]")
        else:
            built_info = (built_info + "\n") + ("#Flanking sequence filter: [Forced OFF]")

        if gapflag ==True:
            built_info = (built_info + "\n") + ("#Min. gap between repeats: [%d bp | aa]"% min_gap)

        if (stdFlag == True) and (dnaFlag == "dna" or dnaFlag == "rna"):
            built_info = (built_info + "\n") + ("#Motif standarization    : [FULL: Complementary + Cyclic]")
        else:
            built_info = (built_info + "\n") + ("#Motif standarization    : [PARTIAL: Cyclic]")

        misa =getMISA(misa_string)                                                      # get the misa format minimum repeat. Note min repeat overwrites misa strings if not zero

        if min_repeat !=0:                                                              # if -min is >0 then overwrite default misa string 
            for x in range(uml+2):
                misa[x] = int(min_repeat)                                               # assing the same minimum repeat for all

        if dict_file == "":
            built_info = (built_info + "\n") + ("#MISA [mininum repeat]   : " +  str(misa))
        else:
            built_info = (built_info + "\n") + ("#Fixed [mininum repeat]  : " +  "["+ str(uml) + "]")

        if (engFlag ==False) and  (dict_file == ""):
            built_info = (built_info + "\n") + ("#Motif dictionary type   : [Construct fast(less accurate) unit motif sequence generator]")
        elif (engFlag ==False) and (dict_file !=""):
            built_info = (built_info + "\n") + ('#Motif dictionary type   : ['+ dict_file + ']')
        else:
            built_info = (built_info + "\n") + ("#Motif dictionary type   : [Construct slow(more reliable) unit motif sequence generator]")
                                  
        print (built_info)                                                              #FINALLY! print the messages combined)

        if (make_primer == True):
            print("Building primers         : [Yes, --> " + primer_saveFileName + "]") 
                                                                 
        if ((fasta_flag == False) and (report_flag==True)):                             # write built_info if teh report atatch is true
            f.write (str(built_info)) 

        # On screen warning messages
        if dict_warning ==True:
        	print("Warning! Dictionary mode is scanning with default values as not all required parameters [-unit, -min and -fix] were supplied!")

        # Pre-run messages  before results real header
        #===============# write the header file into newly created empty output file

        outTXT =("\nsequence_id"+ "\t" + "motif"           + "\t" + "length"             + "\t" + "repeats" + "\t" + "identity(%)" + "\t" + "score" 
                                + "\t" + "seed_missmatch"  + "\t" + "non-seed_missmatch" + "\t" + "start"   + "\t" + "stop"
                                + "\t" + "sequence"        + "\t" + "left"               + "\t" + "right" )
        if (revComp==True):
            outTXT =(outTXT     + "\t" + "ms_revcomp"      + "\t" + "lf_rcomp"           + "\t" + "rf_rcomp")   

        f.write("#Warning: OUTPUT is filtered to contains IMPERFECT repeats only\n") if (align_flag == True and imperfect_only == True) else  None
        f.write("#Warning:The dictionary mode is scanning with default values!")     if (dict_warning == True) else  None 

        f.write(outTXT) if (fasta_flag == False and align_flag == False) else  f.write("\n")                 # don't write header for fasta file output

        print(25 * "__" + "\n")
        #================
        print ("\n#Processing motifs       : May take longer processing time for generating larger (>6) motif lengths")
        
        pat=[]                                                                         # create pattery once before reading the file 
        motif_len_set =[]                                                              # initalize motif length list used
        for pattern_len in range(uFix, uml+1):                                         # pattern_length shoudl start from 0 as python array start with 0          
            print ("\n-|Processing motif length : [%d]"% pattern_len, end="")
            pat.append(re.compile(build_regex(pattern_len, misa[pattern_len-1], imperfect, me_block, engFlag)))  # build and compile re expression len min to max            
            motif_len_set.append(pattern_len) if pattern_len > 0 else None             
        
        print("\n\n#Current Process status : [ Starting sequence scan ]\n")

        for seq_record in SeqIO.parse(handel, "fasta"):
            primer_sum = ("--> " + str(primer_count)) if (primer_count > 0)  else ""      # display cumulatively found primers or ""
            print ("Processing sequence ID  : [", seq_record.id, "]----C.Freq|"+ str(total_ms_count) +"|"+ primer_sum)
            seq=str(seq_record.seq)
            total_seq_scan_length +=len(seq)                                              # get sum of total scanned length 
            total_fasta_seq +=1                                                           # get total numebr of fasta sequence files scanned
            GC_avg  += GC(seq) if (stats_flag == True) and (dnaFlag == "dna" or dnaFlag == "rna") else 0
            results = []
            # make pattern 
            
            tmp_table =[]
            pat_array = 0                                                                 # start reading the pat[] from start 0
            for unit_motif in motif_len_set: #range(uFix, uml+1):   # 1                   # pat[ ] shoud always start with 0 so has to add uFix 
                matches = [[i.start(), i.end(), i.group(), i.group()[:unit_motif]] for i in re.finditer(pat[pat_array], seq)]
                for x in matches:
                    tmp_table.append(x)
                pat_array +=1

            matches=tmp_table   

            for (start, end, micro, seed) in matches:
                    results.append([seq_record.id, start, end, get_std_seed(seed), micro])  # let the motif seed arrange it into most alphabetically  ordered form of all its rotatory forms
            answer= include_imperfect_repeats(results, imperfect, strictFlag, seq, gapflag, min_gap)

            for x in answer:
                ms_type= str(x[3])                                         
                msp=     x[1]                                                        #  have to add +1 during display as python array start from 0
                mep=     x[2]
                mrn=     (int(mep)-int(msp))/int(len(ms_type))                               
                mrn=     round(mrn,3)
                mss=     str(x[4])

                if ( misa[len(ms_type)-1] <= mrn ):                                   # only output result with MSIA features for , match the minimum reepat for the given motif length supplied by MISA
                    #========================= Flanking sequence extract [ok!]

                    if left_ext_len >= msp:
                        LFseq= seq[1:msp]                                             # etract whatever remaning of Left length -1 length (probably at start)
                    else:
                        LFseq= seq[int(msp) -left_ext_len:msp]                             

                    if right_ext_len >= (len(seq)-mep):                               # extract whatever right remaning Length of seq -1 length
                        RFseq= seq[int(mep):len(seq)-1]
                    else:       
                        RFseq= seq[int(mep):int(mep) +right_ext_len]

                    #========================= Calculate score/ mismatches for imperfect repeat
                    scoreList = missmatch_penalty(mss, ms_type, me_block, mmp_val, im_occur_penalty, miss_occur_penalty)
                    score              =  scoreList[0] * (-1)                         # get the final score after processing all, put negative as it is far from perfect, if not zero for perfect
                    non_seed_missmatch =  scoreList[2]        
                    seed_missmatch     =  scoreList[3]       
                    seq_identity       =  scoreList[4]        
                    #========================additional scoring matrix information
                    if ((show_scoreflag == True) and (seq_identity >= sim_value)):    # qualify to show only if results satisfy mrn and fsc parameter after match (makes slow)
                        if (flanking_check(LFseq) == False and flanking_check(RFseq) == False) and (sdout_flag == True):
                            print ("\nComputed score    : ", scoreList[0] * (-1) ) 
                            print ("Base non-match    : ",   scoreList[1] )
                            print ("Non_seed missmatch: ",   scoreList[2] )
                            print ("Seed_missmatch    : ",   scoreList[3] )
                            print ("Seq Identity(%)   : ",   scoreList[4] )
                    #========================= Standarize the seed for both reverse complents
                    org_motif =ms_type                                                                        # conserve original motif type (rotatory equv. ok) for alignment option
                    ms_type = std_motif(ms_type) if ( dnaFlag != "protein" and stdFlag == True) else ms_type  # standarize the seed if condition is true  for -std if not return same ms_type
                    #========================= Write into file and console std-out 
                    outTXT= (     str(seq_record.id) + "\t"  + str(ms_type)                     + "\t" + str(len(ms_type))
                                + "\t" + str(mrn)      + "\t"  + str(seq_identity)              + "\t" + str(score)  
                                + "\t" + str(seed_missmatch)   + "\t" + str(non_seed_missmatch) + "\t" + str(msp+1)   # becuase teh array starts from zero 
                                + "\t" + str(mep)              + "\t" + mss                     + "\t" + LFseq 
                                + "\t" + RFseq )
                    
                    if (revComp==True and (dnaFlag == "dna" or dnaFlag == "rna")):                                      # writes the reverse complementary sequences
                        outTXT =(outTXT + "\t" + seq_revComp(mss) + "\t" + seq_revComp(LFseq)  + "\t" + seq_revComp(RFseq) )

                    # writing  the data to file starts here...with final conditions

                    if ((fsc_flag ==True)  and (seq_identity >= sim_value)):                        # prints the value based on sequence similarity thresthold   and FSC parameter                                                           # check if the flanking sequece contains repeats  for the option (output to the file each new line)     
                        if (flanking_check(LFseq) == False and flanking_check(RFseq) == False):     # as fasta or normal file format
                            f.write ("\n" + outTXT) if (fasta_flag == False and align_flag == False) else None          # write custom file
                            f.write ("\n" +">" + str(seq_record.id) + "|" + str(msp) + "-" + str(mep) + "|"+ str(mrn) + "|" + str(ms_type) \
                                               + "\n" + LFseq + mss + RFseq ) if (fasta_flag == True and align_flag == False) else None                 
                            
                            #===== write alignemnt for only imperfect repeats if true
                            if (imperfect_only == True and score != 0):
                                f.write(get_alignemnt(seq_record.id, mss, mrn, org_motif, msp, mep, score,non_seed_missmatch, seed_missmatch, LFseq, RFseq)) if (align_flag == True and fasta_flag == False) else None
                            elif (imperfect_only == False):
                                f.write(get_alignemnt(seq_record.id, mss, mrn, org_motif, msp, mep, score,non_seed_missmatch, seed_missmatch, LFseq, RFseq)) if (align_flag == True and fasta_flag == False) else None

                            # print on console the results      
                            print   (       outTXT) if (sdout_flag == True) else None           # prints the standard out on console

                            if (make_primer == True):
                                #--------------------make primer using primer3_core If possible using primer 3                              
                                primerID     = seq_record.id + "_"+ ms_type + "_" + str(msp+1)+ "_"+ str(mep)
                                target_len   = str(len(LFseq)) + "," + str(int(len(ms_type) * mrn))            # msp, (repeat length)
                                sequence     = LFseq + mss + RFseq
                                config_data  = primer3_config  (primerID, sequence, target_len, prm_range, prm_opt_size, prm_min_size, prm_max_size, 
                                                                prm_opt_tm, prm_min_tm, prm_max_tm, prm_tm_diff, prm_opt_gc, prm_min_gc, prm_max_gc, prm_poly ,"1")                               
                                written      = write_primer_seq(primer_saveFileName, config_data, sequence, write_header, l_tag, r_tag,ms_type,mrn)
                                write_header = False                         # make it permanenty False next time

                            total_ms_count +=1                          # count total number of mss found
                            #-------------------------------------------
                            if (stats_flag ==True): 
                                total_ms_coverage += len(mss)                                                   # get total coverage length fol all mss
                                motif_freq_byLength [len (ms_type)] +=1
                                motif_coverage_byLength [len (ms_type)] += int(len(ms_type) * mrn)
                                motif_log = (motif_log + "\t" + str(ms_type))
                                imperfect_count += 1 if (imperfect > 0 and score < 0) else 0
                                if (make_primer == True and written != "-1"):
                                    primer_count +=1                          # count primer
                                ms_gc += GC(seq)   if (dnaFlag == "dna" or dnaFlag == "rna") else 0                           # collect GC percenatge of microsatellite
                            
                    elif ((fsc_flag ==False)  and (seq_identity >= sim_value)):                     # prints the value based on sequence similarity thresthold   and FSC parameter                                                                              # writes the output (as fasta or normal) to the file each new line
                        f.write ("\n" + outTXT) if (fasta_flag == False and align_flag == False) else None
                        f.write ("\n" +">" + str(seq_record.id) + "|" + str(msp) + "-" + str(mep) + "|"+ str(mrn) + "|" + str(ms_type) \
                                                + "\n" + LFseq + mss + RFseq ) if (fasta_flag == True and align_flag == False) else None
                        #===== write alignemnt for only imperfect repeats if true
                        if (imperfect_only == True and score != 0):
                            f.write(get_alignemnt(seq_record.id, mss, mrn, org_motif, msp, mep, score,non_seed_missmatch, seed_missmatch, LFseq, RFseq)) if (align_flag == True and fasta_flag == False) else None
                        elif (imperfect_only == False):
                            f.write(get_alignemnt(seq_record.id, mss, mrn, org_motif, msp, mep, score,non_seed_missmatch, seed_missmatch, LFseq, RFseq)) if (align_flag == True and fasta_flag == False) else None
                        #===== print output on console
                        print   (       outTXT) if (sdout_flag == True) else primer_count           # prints the standard out on console

                        if (make_primer == True):
                            #--------------------make primer using primer3_core If possible using primer 3
                            primerID     = ms_type+"_"+str(msp)+ "_"+ str(mep)
                            target_len   = str(len(LFseq)) + "," + str(int(len(ms_type) * mrn))            # msp, (repeat length)
                            sequence     = LFseq + mss + RFseq
                            config_data  = primer3_config  (primerID, sequence, target_len, prm_range, prm_opt_size, prm_min_size, prm_max_size, 
                                                            prm_opt_tm, prm_min_tm, prm_max_tm, prm_tm_diff, prm_opt_gc, prm_min_gc, prm_max_gc, prm_poly ,"1")
                            written      = write_primer_seq(primer_saveFileName, config_data, sequence, write_header, l_tag, r_tag, ms_type, mrn)
                            write_header = False
                            
                        total_ms_count +=1                                                           # count total number of mss found
                            #--------------------------------------------
                        if (stats_flag ==True):                                                             
                            total_ms_coverage += len(mss)                                            # get total coverage length fol all mss
                            motif_freq_byLength [len (ms_type)] +=1
                            motif_coverage_byLength [len (ms_type)] += int(len(ms_type) * mrn)          
                            motif_log = (motif_log + "\t" + str(ms_type))                            # only claculate if flag is true
                            imperfect_count += 1 if (imperfect > 0 and score < 0) else 0
                            if (make_primer == True and written != "-1"):                            # count primers
                                primer_count +=1
                            ms_gc += GC(mss) if (dnaFlag == "dna" or dnaFlag == "rna") else 0                      # collect GC eprcentage of mss

    #--------------------------------------------------------------------------------------------------------------------------
    print ("\n----------[end of scan]----------")

    if (stats_flag == True):                                      # printout the report on stdout 
        lst =[]
        report_data =""                                           # holds teh final report information
        for x in get_report (motif_log,dnaFlag):
            lst.append((x[0],x[1]))
        lst=sorted(lst)
        report_data += "\n" + ("Report: Basic statistics")
        report_data += "\n" + ("[Motif]" + 2 * "\t" + "[Frequency]")
        report_data += "\n" + ("--------------------------------------")
        lst.sort(key=lambda x: x[1])
        for x in lst:
            report_data += "\n" + (x[0]  + 2 * "\t" + str(x[1]))
        report_data += "\n" +     ("--------------------------------------")
        report_data += "\n" +     ("Sequence length scanned  (bp|aa):\t%d"% total_seq_scan_length)
        report_data += "\n" +     ("Number of sequence file(s)      :\t%d"% total_fasta_seq)

        if (dnaFlag == "dna" or dnaFlag == "rna"):
            report_data += "\n"     + ("Average G+C Percentage          :\t"+ str(round(GC_avg/total_fasta_seq, 4)) )
            if total_ms_count > 0:
                report_data += "\n" + ("Average G+C Percentage in repeat:\t"+ str(round(ms_gc/total_ms_count  , 4)) )

        report_data += "\n" +     ("Total repeat motif(s) identified:\t%d"% total_ms_count)
        if (imperfect > 0):
            report_data += "\n" + ("Imperfect repeats identified    :\t%d"% imperfect_count)
        report_data += "\n" +     ("Relative density       (No./Mbp):\t%f"% round(1000000 * total_ms_count/total_seq_scan_length,4))
        report_data += "\n" +     ("Repeat motif(s) coverage (bp|aa):\t%d"% total_ms_coverage)
        report_data += "\n" +     ("Percentage coverage by repeats  :\t%f"% round(100 * total_ms_coverage/total_seq_scan_length,4))  
        for x in range(uFix, uml+1):
            report_data += "\n" + ("Freq., coverage(bp), density (No./Mbp) of motif[" + str(x) +"]:\t" + \
                            str(motif_freq_byLength[x]) + "\t"+ str(motif_coverage_byLength[x]) + "\t"   + \
                            str(round(1000000 * motif_coverage_byLength[x]/total_seq_scan_length,4)) )  # display in mss/million base pair
        report_data += "\n" + ("--------------------------------------") 

        print(report_data)                                                                  # print report
        print ("Total primer pair(s) generated  :\t", primer_count) if (make_primer == True ) else None

    
    if ((fasta_flag == False) and (report_flag==True)):                                     # open the report file and append the report at end
        with open(outfile,"a+") as f:
            f.write ("\n----------[end of scan]----------") 
            f.write (str("\n" + report_data) ) 

    print ("\nThank you for using OSTRFPD...best wishes\n" )
    print ("\nOperational time elapsed:", abs(round(start_time - time.time(),3)), "s")      # get elapsed time for scan
    
    return lst                                                                               # well... alldone here.... end of the function and return

#=====================Declaring command line arguments type
parse = argparse.ArgumentParser()

parse.add_argument("-input" , type = str,  help = "Full pathname of input file (i.e. source FASTA file). If not supplied, the [input.fasta]\
                                                   file will be searched by default. Single, multi-FASTA and gunzip-compressed (.gz) FASTA are supported for direct scan without unzipping")

parse.add_argument("-output", type = str,  help = "Full pathname of output result file. If not supplied, the < source filename + _res.txt > will be used as default output report file. \
                                                   The default output file is tab-delimited plain text file containing multiple parameters (e.g.: unit motif, repeat number, start stop \
                                                   position, flanking sequences). The output file can also be configured to contain built information and general motif statistics and categorization")  

parse.add_argument("-unitmax"  , default = 6, type = int, help = "Input type: positive integer. Sets the maximum unit motif length to be searched.\
                                                    \nIf used with ‘-min’ argument only MIN single unit fixed will be scanned.\
                                                    \nIn Fast (default) search mode, the maximum values for DNA,RNA and Amino acid unit motifs \
                                                      are 10, 10 and 3 bp|aa respectively.\
                                                    \nIn Accurate search mode, the maximum values for DNA, RNA and Amino acid unit motifs are 6, 6 and 3 bp|aa respectively.\
                                                    \nIn Dictionary-based search mode, the maximum values for DNA, RNA and Amino acid unit motifs is 30 bp|aa.")

parse.add_argument("-unitmin"  , default = 1, type = int, help = "Input type: positive integer. This option sets the minimum unit motif length \
                                                     to be searched. Default values is 1. The UNITMIN can range from 1 to UNITMAX.If UNITMIN \
                                                     and UNITMAX are same, then single fixed length unit motif will be search. Use of [-fix true] \
                                                     argument will overwrite the ‘-unitmin UNITMIN’ setting and forcefully set the UNITMIN = UNITMAX (i.e. single fixed unit motif.")

parse.add_argument("-min"   , default = None, type= int
                            , help =               "Input type: Positive integer. Sets the minimum repeats (copy number) threshold for selection of \
                                                    tandemly repeated sequences. This option will overwrite MISA-formatted minimum repeat settings\
                                                     (if present). By default the MIN value is 0 (i.e. OFF).")

parse.add_argument("-imperfect", type =int, help=  "Maximum distance (aa|bp) within which the mismatch (indel) containing tandem repeats but with same\
                                                    unit motif (seed) will be combined as single imperfect repeat. The value for argument ‘–imperfect’ > 0,\
                                                    automatically turns ON search for imperfect repeats as associated default score values [i.e. –imop, -mop, -mmp]")

parse.add_argument("-lflank", type =int, default =85,   help = "Input type: Positive integer. Extracts the left flanking sequence. This option auto-truncates the \
                                                   flanking sequence region if the left side sequence is smaller than supplied length (bp | aa) of flanking sequence. Default value is 85.")

parse.add_argument("-rflank", type =int, default =85,   help = "Input type: Positive integer. Extracts the right flanking sequence. This option auto-truncates the \
                                                   flanking sequence region if the right side sequence is smaller than supplied length (bp | aa) of flanking sequence. Default value is 85.") 

parse.add_argument("-exclude", default =          "false", choices =[None, "true", "false" ] 
                             , help =             "Restricted or Strict mode. Ignores non-seed (aa|bp) containing imperfect microsatellites even if the non-seed mismatch occurred \
                                                   within pre-declared imperfection range. Default value is None or false.")

parse.add_argument("-imop"   , type =int   
                             , help =              "Input type: positive integer. Penalty score for the first non-seed mismatch occurrence (added only once for an\
                                                    imperfect microsatellite independently on top of ‘-mop’ and ‘-mss’ parameter values). Default value is 3 and \
                                                    only applied if value for ‘-imperfect’ > 0 (i.e. IMPERFECT mode is set ON).")

parse.add_argument("-mop"   , type =int   
                            ,               help ="Input type: positive integer. Penalty score which is applied only once for the first (starting mismatch) occurrence\
                                                   to a given imperfect microsatellite. Default value is 2 only applied if value for ‘-imperfect’ > 0 (i.e. IMPERFECT mode is set ON).")

parse.add_argument("-scan"  , default = "dna" 
                            , choices =[None, "dna","rna" , "protein" ]
                            ,               help ="Sets the expected source sequence type to DNA, RNA, or PROTEIN. Default value is set to dna") 

parse.add_argument("-mmp"   , type =int   
                            , help=               "Input type: positive integer. Set the penalty for each mismatch (seed or non-seed type).\
                                                   It is added on top of mismatch starting penalties supplied by '-impo' and/or '-mop' if present. \
                                                   Default value is 2 and only applied if value for ‘-imperfect’ > 0 (i.e. IMPERFECT mode is set ON).")

parse.add_argument("-rcomp" , default = None, choices =[None, "true","false"]
                            , help =  "Attaches the reverse complement sequence for DNA or RNA repeats including flanking sequence. Default value is None or false.")

parse.add_argument("-misa"  , default = None, type =str
                            , help = "MISA-formatted number series to input different minimum repeat number for different motif length. \
                                     [e.g: -misa 14,7,5,4,4 for minimum repeat number of microsatellites with unit motif length of 1,2,3,4 and 5, respectively].\
                                     Note that the minimum repeat values for misa should be supplied as [1,2,3,...,UNITMAX-1,UNITMAX] even\
                                     if minimum unit motif length (UNITMIN)>1. The misa is automatically set as NULL if [-min] > 0. \
                                     Default value of misa string for DNA or RNA is 14,7,5,4,4,4 and Protein is 7,5,3,3 , respectively.\
                                     Whenever possible do avoid using ‘-misa’ and ‘-min’ arguments together to make syntax clear for operation.") 

parse.add_argument("-eng"   , default = "false", choices =["true","false"]
                            , help =  "Sets the unit motif sequence generator engine type for accurate fast or accurate search. Default value is [false]\
                                       for fast but less efficient motif pattern search. If set ‘true’ then uses slow yet accurate engine.")

parse.add_argument("-sshow" , default = "false", choices =["true","false"]
                            , help =  "Sets the display output ON/OFF for the scoring matrix on console screen during scan for each tandem repeats \
                                       identified. Has no effect on output result file. Default value is false. This option is forced turn OFF if [-sdout] is false") 

parse.add_argument("-fsc"   , default = None, type =str
                            , help = "Used as positive integers in the format ‘m,n’ (Typical use: [-fsc 1,5] or [-fsc 1,6]). Default status is OFF or None.\
                                      Filters out the microsatellites if a minimum of 'm' unit motif length tandem repeat of minimum 'n' repeat number (copy number) \
                                      is found in Left or Right flanking sequence. Note: The option is useful for removing microsatellites from results \
                                      whose low-numbered repeats may cause problem in primer design.")

parse.add_argument("-pfname", default = None, type =str
                            , help = "Full pathname of primer result output file to be saved. If not supplied, by default the output name will be used \
                                      as [source filename + ‘_prm.txt’].The output is a tad-delimited plain text file. The primer output file will \
                                      contain Sequence IDs of primer derived from the original sequences from which the primers were designed. \
                                      In addition to left (forward) and right (reverse) primers, the corresponding standardized unit motifs, \
                                      motif copy number, Tm and complete sequence (Left flanking+ microsatellite + right flanking) will be tabulated in the primer result file.") 

parse.add_argument("-fix"   , default ="false", choices =["true","false"]
                            , help = "Sets a single fixed minimum repeats (copy number) for all unit motifs. This option will overwrite MISA settings (if present). Default value is false.") 

parse.add_argument("-primer", default ="false", choices =["true","false"]
                            , help = "Sets the flag to make primer using Primer3. Default value is false.") 

parse.add_argument("-std"   , default ="true", choices =[None, "true","false"]
                            , help = "Standardizes the unit motif (seed) name for display and report in either \
                                      Partial:[cyclic equivalent motifs] or Full: [reverse complement + cyclic equivalents ].\
                                      Default is set to Full standardization. See citation paper for details. \
                                      General statistics results for motif characterization will also be formatted based on the \
                                      type of motif standardization. For Protein sequence, as there are no complement strands, \
                                      only partial standardization is utilized for amino acid unit motif characterization. ") 

parse.add_argument("-dict"  , default ="false", type=str  
                            , help = "[Usage: -dict filename]. Dictionary file is a plain text file contaning list of unit motifs of\
                                       same type molecules (DNA, RNA, or amino acids) each separated by new line. Dictionary file is \
                                       used to exclusively search custom user-supplied unit motif sequences. If used, only the unit \
                                       motifs listed in the dictionary file (e.g. dict_dna.txt) will be searched for a single fixed \
                                       (default: -fix true) length (default: -unit 8) unit motifs of maximum upto 30 bp and fixed \
                                       minimum repeat number (default: -min 4). This mode [-dict] is expected to be used together \
                                       with [-unit, -fix and -min] parameters. Please view details in software documentation manuel.") 

parse.add_argument("-fasta" , default ="false", choices =["true","false"]
                            , help = "Sets the output result to FASTA formatted file. File description, motif statistics or built\
                                      headers is NOT available for this Mode. If output filename is not supplied source filename + ‘_res.txt’ \
                                      will we written as default output. Default values if false.") 

parse.add_argument("-sdout" , default ="false", choices =["true","false"]
                            , help = "Sets the standard result display flag for console or terminal. If true, displays each identified \
                                      repeat result on console or terminal screen. If set false, only the Sequence IDs being search along \
                                      with cumulative frequency (C.Freq |total|) of identified repeats and primer designed (Total primers designed) \
                                      will be shown. Has no effect on output result file. If False, Default value is false.")   

parse.add_argument("-gap"   , default = 0, type= int
                            , help = "Input type: Positive integer.  Sets minimum gaps ( aa | bp ) between two consecutive repeat motifs.\
                                      Please note that all consecutive repeats identified within the range will be discarded even if it \
                                      belongs to same motif type. Default value is zero. The default value is 0.")

parse.add_argument("-sim"   , default = 50, type= int
                            , help = "Input type: positive integer. Sets the minimum similarity threshold value (in percentage, e.g: -sim 50 for 50 percent)\
                                      for the results to be accepted on top of other selection criteria supplied. Similarity basically refers to percentage\
                                      of complete unit motifs present in the imperfect repeat sequence compared to its near equivalent perfect repeat.\
                                      All perfect repeats have similarity value of 100 percent (See citation paper for details). Default value is 50.")

parse.add_argument("-stats" , default ="true", choices =["true","false"]
                            , help = "Displays and append basic statistics of the results to output file \
                                      (option not applicable for FASTA output)")

parse.add_argument("-report", default ="true", choices =["true","false"]
                            , help = "Appends Built header on top and Basic statistics report at the end of the output file (not applicable for FASTA output)")

parse.add_argument("-align" , default ="false", choices = [None, "true","false"]                                       
                            , help = "Output the report file in alignment form (not applicable for FASTA or Custom format)")

parse.add_argument("-autoexit", default ="true", choices =["true","false"]
                              , help = "Sets auto exit flag. If true, auto exits after completion of task or pauses the console \
                                      screen if false (applicable mainly for windows system GUI launcher). Default value is false.") 

parse.add_argument("-ltag"  , default = None, type =str
                            , help = "Attaches the left tag (e.g: 6FAM-) for the left primer generated by Primer3.\
                                       Typical use (‘-ltag 6FAM-’). Default value is None.") 

parse.add_argument("-rtag"  , default = None, type =str
                            , help = "Attach the right tag (e.g: gtgtctt-) for the right primer generated by Primer3.\
                                      Typical use [-rtag gtgtctt- ]. Default value is None") 

parse.add_argument("-imalign", default ="false", choices = [None,"true","false"] 
                             , help = "Make the alignment output file only contain imperfect alignment results after satisfying \
                                       all other applicable conditions. The alignment option will be automatically activated.") 

parse.add_argument("-gui"   , default ="false", choices = [None,"true","false"] 
                            ,  help = "Open the OSTFRPD in GUI interface to input the configuration parameters.\
                                       Using this option will OVERWRITE all other CLI arguments supplied and launches OSTRFPD \
                                       is fresh default GUI mode. Default value is false.")

parse.add_argument("-prng"  , default = "150-300", type =str
                            , help = "Input type: String. Used to set min-max (e.g: 150-300), the minimum and maximum range\
                                      (in bp) of the output primers product (amplicon size). Default value is 150-300.") 

parse.add_argument("-posz"  , default = "20", type =str
                            , help = "Input type: positive integer. Used to set the optimum primer length (in Bp). Default value is 20")

parse.add_argument("-pmisz" , default = "17", type =str
                            , help =  "Input type: positive integer. Used to set the minimum primer length (in Bp). Default value is 17") 

parse.add_argument("-pmxsz"  , default = "26", type =str
                            , help = "Input type: positive integer. Used to set the maximum primer length (in Bp). Default value is 26") 

parse.add_argument("-potm"  , default = "60", type =str
                            , help = "Input type: positive integer. Used to set the optimum primer Tm (C). Default value is 60") 

parse.add_argument("-pmitm"  , default = "58", type =str
                            , help = "Input type: positive integer. Used to set the minimum primer Tm (C). Default value is 58") 

parse.add_argument("-pmxtm" , default = "63", type =str
                            , help = "Input type: positive integer. Used to set the maximum primer Tm (C). Default value is 63") 

parse.add_argument("-ptmdiff"  , default = "6", type =str
                            , help = "Input type: positive integer. Used to set the maximum primer Tm (C) difference. Default value is 6")

parse.add_argument("-pogc"  , default = "55", type =str
                            , help = "Input type: positive integer. Used to set the optimum GC content (given in percentage). Default value is 55") 

parse.add_argument("-pmigc" , default = "20", type =str
                            , help = "Input type: positive integer. Used to set the minimum GC content (given in percentage). Default value is 20") 

parse.add_argument("-pmxgc" , default = "80", type =str
                            , help = "Input type: positive integer. Used to set the maximum GC content. Default value is 80") 

parse.add_argument("-pmpoly", default = "3", type =str
                            , help = "Input type: positive integer. Used to set the maximum Poly-X's in primer. Default value is 3") 

parse.add_argument("-v", default = "", type =str
                            , help = "(O)mni (S)hort (T)andem (R)epeat (F)inder & (P)rimer (D)esigner: OSTRFPD Version 0.01. Developed by Mathema VB and Imwong M - 2018.(under review...)" ) 

args = parse.parse_args()        # set the argument parsing setting to class variable args

#=======================================
if (args.primer == "true"):             #must be declared before the pfname parameter           
    make_primer = True                                                 
else:
    make_primer =  False
#========================Argument process
if (args.input== None):
    targetfile = "input.fasta"
    outfile    = "output.txt" 
else:
    targetfile = str(args.input)
    outfile    = (str(args.input) + "_res.txt") 


if (args.output != None):
    outfile = str(args.output) + "_res.txt" 

elif (args.output== None and args.input != None):
    outfile = (str(args.input) + "_res.txt") 

elif (args.output== None and args.input== None):
    outfile    = "output.txt" 


if args.pfname != None:
    primer_saveFileName = str(args.pfname)

elif (args.input != None):
    primer_saveFileName = (str(args.input) +"_prm.txt") 

else:
    primer_saveFileName = "primer_results.txt"


if (args.unitmax == None):                                               
    unit_motif_len= 3                                                 # default is three (for fast)
else:
    unit_motif_len= int(args.unitmax)

if (args.unitmin == None):                                           
    unit_motif_min_len= 1                                             # default is three (for fast)
else:
    unit_motif_min_len= int(args.unitmin)

if (args.min == None or args.min =="0"):
    min_repeat= 0
else:
    min_repeat= int(args.min)                                          # Because: Total Reepat No = total occurance - unit occurance
    misa_string =""                                                    # overwrite misa string

if (args.imperfect == None):
    imperfect_repeate = 0                                              # change later
else: 
    imperfect_repeate= int(args.imperfect)

if (args.lflank == None):
    left_ext_len =  85                                                # left flanking sequence to extract (default 85)
else:
    left_ext_len = int(args.lflank)

if (args.rflank == None):
    right_ext_len = 85                                                 # right flanking sequence to extract (default 85)
else:
    right_ext_len = int(args.rflank)

if (args.mmp == None):
    mmp_val = 1                                                        # default right now is set to one
else:
    mmp_val = int(args.mmp)

if (args.scan == None or args.scan == "dna"):                         # select among DNA/ RNA or protein (default DNA)                              
    dnaFlag = "dna"
elif (args.scan == "rna"):
    dnaFlag = "rna"
elif (args.scan == "protein"): 
    dnaFlag = "protein"
else:
    dnaFlag = "dna"

if (args.gap == None):                                  
    gapflag =False
    min_gap =0
else:
    gapflag =True
    min_gap =int(args.gap) 

if (args.sdout == None):                                              # settings for standard display out for results                           
    sdout_flag = True
elif (args.sdout == "false"):
    sdout_flag = False
else:
    sdout_flag = True

if (args.stats == None or args.stats == "true" ):                                              # calculate and show statistics                           
    stats_flag = True
elif (args.stats == "false"):
    stats_flag = False
else:
    stats_flag = True

if (args.exclude == None or args.exclude == "false"):              
    strictFlag = False
else:
    strictFlag = True

if (args.std == None or args.std == "true" ):                         # standarize the seed for its reverse complement            
    stdFlag = True
elif args.std == "false":
    stdFlag = False
else:
    stdFlag = True

if (args.align == None or args.align == "false"):                     # set the alignment output file style type      
    align_flag = False
elif args.align =="true":
    align_flag = True
else:
    align_flag = False

if (args.report == None or args.report == "true" ):                   # condition to attach report file       
    report_flag = True
elif (args.report == "false"):
    report_flag = False
else:
    report_flag = True

if (args.fasta == None):                                              # fasta file flag settings
    fasta_flag = False
elif (args.fasta == "true"):
    fasta_flag = True
else:
    fasta_flag = False

if (args.imop == None):
    im_occur_penalty =   3                                            # defualt is right now set to 10
else:
    im_occur_penalty = int(args.imop)

if (args.mop == None): 
    miss_occur_penalty = 2                                           # default set to 10 [ is auto zero if im_occur_penalty=0]
else:
    miss_occur_penalty = int(args.mop)

if (args.rcomp == None or args.rcomp == "false"):                    # reverse complement of the entire sequence
    revComp = False                                                  
elif (args.rcomp == "true" ):
    revComp = True

if (args.fix == "true"):                                             # set the fix number (min =1 as pat[] beings fro m[1]) of the unit to scan, will give error if not used unit
    uFix_flag = True                                                
elif (args.fix =="false" or args.fix == None):
    uFix_flag = False

if (args.misa == None or args.misa.strip() ==""):                                              # process the MISA string
    misa_string = ""                                                   
else:
    misa_string = str(args.misa)

if (args.min != None):
    misa_string = "" 

if (args.imalign == None or args.imalign == "false"):                # process the imperfect micorsatellite only in output
    imperfect_only = False                                                 
else:
    imperfect_only = True
    align_flag = True                                                 # automatically make fasta file output true                   

if (args.eng == None):                                                # process teh MISA string
    engFlag  = False                                                  # use fast but less reliable engine by default                                       
elif (args.eng == "true"):                                            # use more reliable engine
    engFlag = True                                            
elif (args.eng == "false"):                                           # use the less reliabel engine by settings
    engFlag = False 
else:
    engFlag = False                                                   # just in case of misstye use fast engine

if (args.sshow == None or args.sshow == "false" ):                    # process teh MISA string
    show_scoreflag  = False                                                
elif (args.sshow == "true"):
    show_scoreflag  = True

if (args.fsc == None or args.fsc == 'none' ):                                                 # process teh MISA string
    fsc_flag = False
    fsc_string = ""                                                   # e.g: 1;6 which sets minimum unit motif of 1 and minimum 6 repeats (Note: must need to set it in global varaible even if False to put in recursive fucntion)
else:                                                                 # provide string in format x;y where x is minimum motif and y is m
    fsc_flag =   True
    fsc_string = str(args.fsc)                                        # set the fsc parameetr if the -fsc  is supplied 

if (args.sim != None):                                                # sets the similarity value thresthold                
    sim_value =abs(int(args.sim))
elif (args.sim == None or args.sim == False):
    sim_value =100
else:
    sim_value =100

if (args.ltag != None):                                             # for left priemr tag argument
    l_tag = str(args.ltag) 
else:
    l_tag =""

if (args.rtag != None):                                             # for right priemr tag argument
    r_tag = str(args.rtag)
else:
    r_tag =""

if (args.gui == None or args.gui == "false" ):                       # for GUI interface launch
    gui_flag  = False
else:
	gui_flag  = True                                                 # just in case 

dict_warning = False                                                  # initalize dictionary warning
if (args.dict == None or args.dict == "false"):                       # DICTIONARY FILE scan mode
    dict_file = ""
    dict_warning =False
else:
    dict_file =str(args.dict)
    if ((args.fix == "false" or args.fix == None) or (args.min == 0 or args.min == None)):  # Set default dictionary parameters
    	uFix_flag = True
    	min_repeat= 4
    	unit_motif_len = 8
    	dict_warning =True                                               
    
if (args.prng != None):                                             # for range
    prm_range = str(args.prng)
else:
    prm_range ="150-300"

if (args.posz != None):                                             # for range
    prm_opt_size = str(args.posz)
else:
    prm_opt_size ="20"

if (args.pmisz != None):                                             # for range
    prm_min_size = str(args.pmisz)
else:
    prm_min_size ="17"

if (args.pmxsz != None):                                             # for range
    prm_max_size = str(args.pmxsz)
else:
    prm_max_size ="26"

if (args.potm != None):                                             # for range
    prm_opt_tm = str(args.potm)
else:
    prm_opt_tm ="60"

if (args.pmitm != None):                                             # for range
    prm_min_tm = str(args.pmitm)
else:
    prm_min_tm ="58"

if (args.pmxtm != None):                                             # for range
    prm_max_tm = str(args.pmxtm)
else:
    prm_max_tm ="63"

if (args.ptmdiff != None):                                             # for max tm diff
    prm_tm_diff = str(args.ptmdiff)
else:
    prm_tm_diff ="6"

if (args.pogc != None):                                             # for range
    prm_opt_gc = str(args.pogc)
else:
    prm_opt_gc ="55"

if (args.pmigc != None):                                             # for range
    prm_min_gc = str(args.pmigc)
else:
    prm_min_gc ="20"

if (args.pmxgc != None):                                             # for range
    prm_max_gc = str(args.pmxgc)
else:
    prm_max_gc ="80"

if (args.pmpoly != None):                                             # for range
    prm_poly = str(args.pmpoly)
else:
    prm_poly ="3"

#-----------------------------Primer flanking length
if (make_primer == True and (left_ext_len ==0 and right_ext_len == 0) ):
    left_ext_len  =85
    right_ext_len =85

if (args.v == "V" or args.v == "ver"):
	print ("\nOmni Short Tandem Repeat Finder & Primer Designer (OSTRFPD)"  +
		   "\n\nVersion 0.01. Developed by (under review) - 2018" +
		   "\nDepartment of Tropical Medicine, Mahidol University, Bangkok, Thailand" +
		   "\ncontact : Mathema VB, Dondorp AM, Imwong M (2019)" +
           "\n[ use '-gui true' agrument for starting the OSTRFPD in GUI mode ]" +
		   "\n[ use '-h' or '--help' argument to view detailed command arguments ] ")
	sys.exit()    

''' UNCOMMENT THE MODULE for using  motif_chart(lst) Starting here
#====================== OPTIONAL! Code for GENERATING MOTIF GRAPH CHART for visualization
def motif_chart(lst):           # Optional module: Plot bar chart for motif destribution using matplot_lib 
                                # Please use compatible version of Matplot (ver. 2.2.2) or numpy (ver 1.14.5) libs if error.
    print( "Search complete. Generating bar chart from Motif destribution data")
    import matplotlib.pyplot as plt
    import math, numpy as np 

    fig, ax = plt.subplots()
    ax.set_ylabel('Repeat frequency (log10)' )
    ax.set_xlabel('motif sequence')
    ax.set_title ('Destribution of tandem repeats')
    xs, ys = [*zip(*lst)]
    plt.bar(xs, ys)
    
    ax.set_yscale('log')        # display frequency in log 10 as teh difefrences can be huge
    if len(xs) > 10:            # rotates the text x-label to 90 if more than 10 motif category to display   
        for tick in ax.get_xticklabels():
            tick.set_rotation(90) 
    plt.show()
#======================
Uncomment for using  motif_chart(lst) Starting here uptill 
''' 


#-------------------------------------------------------------------MAIN PROGRAM LAUNCHES HERE      
if __name__ == "__main__":

    #gui_flag  = True                                               # UNCOMMENT  to start in GUI mode all time 

    # GUI mode running code block
    if  gui_flag  == True:
        app = QtWidgets.QApplication(sys.argv)
        window = MyApp()  
        window.show()
        warnings.filterwarnings("ignore")                           # supress unnecessary/minor warnings in console output  (version, depreciation stuffs like that)
    #==========Main program control begins here
    print("\nPython version:\n", sys.version)
    print("\n"+  
    "#########  ######  ######## ########  ######## ########  ########  "+ "\n" +   
    "##     ## ##    ##    ##    ##     ## ##       ##     ## ##     ## "+ "\n" +    
    "##     ## ##          ##    ##     ## ##       ##     ## ##     ## "+ "\n" +     
    "##     ##  ######     ##    ########  ######   ########  ##     ## "+ "\n" +     
    "##     ##       ##    ##    ##   ##   ##       ##        ##     ## "+ "\n" +    
    "##     ## ##    ##    ##    ##    ##  ##       ##        ##     ## "+ "\n" +     
    "######### #######     ##    ##     ## ##       ##        ########  "+ "\n") 

    print (".-----------------------------------------------------------------.")
    print ("| (O)MNI (S)HORT (T)ANDEM (R)EPEAT (F)INDER & (P)RIMER (D)ESIGNER |")
    print ("|        <__Mathema VB, Dondorp AM, Imwong M (2019)__>            |")
    print ("|Center for Tropical disease, Mahidol University,Bangkok Thailand |")
    print ("|_________________________________________________________________|")
    print ("\n")
    #===========================================
    write_header = True                                           # to write the header for first occurance


    #-------------------------------------------------------------#Command line interface (CLI) mode code block
    if  gui_flag  == False:
        warnings.filterwarnings("ignore")                         # ignore version-rleated or unnecessary minor warnings

        if True: #(args.input == None):
            retval =Search_Microsatellite (targetfile, outfile, unit_motif_len, unit_motif_min_len, min_repeat, uFix_flag, imperfect_repeate
                                                , left_ext_len, right_ext_len, im_occur_penalty, miss_occur_penalty 
                                                , misa_string, mmp_val, strictFlag, dnaFlag, engFlag, stats_flag 
                                                , show_scoreflag, fsc_flag, fsc_string, fasta_flag, sdout_flag, gapflag
                                                , min_gap, revComp, report_flag, align_flag, sim_value
                                                , dict_file, stdFlag, write_header, make_primer, l_tag, r_tag, imperfect_only
                                                , prm_range, prm_opt_size, prm_min_size, prm_max_size, prm_opt_tm, prm_min_tm
                                                , prm_max_tm,prm_tm_diff, prm_opt_gc, prm_min_gc, prm_max_gc, prm_poly)

    if  (args.autoexit == "false"):                                   # Pause the output screen on processing end ( for Windows if controleld by other app)
        os.system('pause') if os.name == "nt" else None               # Pause the consol from auto-exiting in windows
    sys.exit(app.exec_()) if gui_flag == True else None               # only if the GUI is true