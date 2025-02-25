from engine_ui import Ui_MainWindow
from PyQt5 import QtWidgets
from BioSeq import *

import sys

class DNAENGINE:
    def __init__(self, *args, **kwargs):
        """Creating accessible window we can 'talk' to"""
        super().__init__(*args, **kwargs)
        self.app = QtWidgets.QApplication(sys.argv) # QApplication constructor accessing arugments
        self.MainWindow = QtWidgets.QMainWindow() # Window constructor

    def setup(self):
        """Load Ui file that we generated"""
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.MainWindow)
        self.ui.pushButton.clicked.connect(self.DNA_analysis)

    def DNA_analysis(self):
        """
        run and display everything from Bio_Seq functions after processing.
        \n - need codon frequency option
        """
        input_sequence = self.ui.textEdit.toPlainText().replace("\n", "").strip().upper()

        if not input_sequence:
            self.ui.textBrowser.setText("No sequence entered")
            return

        seq = BioSeq(input_sequence)
        if not seq.validate(): # some error handling in case user puts in something non biotype related like "error_test"
            error_message = (f"{"=" * 60} \n"
                    "Sequence is Invalid\n"
                    "Check for anything other than these characters -> [ATCG]\n"
                    f"{"=" * 60} \n")
            self.ui.textBrowser.setText(error_message)
        else:
            info = (
            f"{"=" * 60} \n"
            f"{seq.get_seq_info()}\n"
            f"{"=" * 60} \n"
            f"[Nucleotide Frequency]: {seq.nuc_freq()}\n"
            f"[RNA Transcribe]: {seq.transcription()}\n"
            f"[Reverse Complement]: {seq.reverse_complement()}\n"
            f"{"=" * 60} \n"
            f"[GC Content]: {seq.gc_content()}\n"
            f"[GC Subsection Content (k=20)]: {seq.gc_subsec()}\n"
            f"{"=" * 60} \n"
            f"[AA Chain]:\n{seq.translate_seq()}\n"
            f"[All Open Reading Frames]:\n{seq.gen_reading_frames()}\n"
            f"{"=" * 60} \n"
            f"[All Possible Proteins]:\n{"\n".join(seq.all_orf_proteins())}\n"
            f"{"=" * 60} \n"
                )

            self.ui.textBrowser.setText(info)


    def run(self):
        """Event loop. Listens to inputs and whatever comes from OS"""
        sys.exit(self.app.exec_())

    def display(self):
        """Renders window on screen"""
        self.MainWindow.show()



