import sys
import os
from graphic_elements.desplegable_menu_1 import DropdownWindow1
from graphic_elements.desplegable_menu_2 import DropdownWindow2
from graphic_elements.desplegable_menu_3 import DropdownWindow3
from graphic_elements.desplegable_menu_4 import DropdownWindow4
from graphic_elements.label_text import CustomLabel
from graphic_elements.py_toggle import CustomCheckBox
from graphic_elements.up_button_triangle import CustomTriangleButton
from graphic_elements.label_for_path import CustomPathLabel
from graphic_elements.MyCombo import CustomCombo
from graphic_elements.up_buttons import CustomButton
from graphic_elements.label_titles import CustomTitleLabel
from graphic_elements.Line_edit import CustomLineEdit
from graphic_elements.Dialog_for_receptors import CheckableOptionsDialog
from codes import Create_pdbqt_by_smi_meeko
from codes import Create_pdbqt_by_txt_meeko
from codes import Create_pdbqt_by_sdf_meeko
from codes import Create_pdbqt_by_txt_obabel
from codes import Create_pdbqt_by_sdf_obabel
from codes import Create_pdbqt_by_smi_obabel
from codes import Create_protein_pdbqt_by_obabel
from codes import Docking_with_smina
from codes import Redocking_with_smina
from codes import Convert_pdbqt_to_mol2
from codes import Calculate_RMSD
from codes import Best_screening
from codes import Consensus_model_for_smi
from codes import Consensus_model_for_txt
from PySide6.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QFrame, QHBoxLayout, QStackedWidget, QPushButton, QSpacerItem, QSizePolicy, QProgressBar, QFileDialog, QMessageBox
from PySide6.QtGui import QIcon, QPixmap, QRegularExpressionValidator, QIntValidator, QDoubleValidator
from PySide6.QtCore import Qt, QRegularExpression


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # "------------------------- ESTABLISH THE MAIN ATTRIBUTES -------------------------"


        # Simulated clicks between windows
        self.ws_1 = False
        self.ws_2 = False
        self.ws_3 = False
        self.ws_4 = False


        # ------------------------- MAIN WINDOW SETTINGS -------------------------"

        
        self.setWindowTitle("DockPipeGui v. 0.1.14.") # Establish the Main Window title
        self.resize(1400, 800)  # Establish main window starting size


        # ------------------------- MAIN FRAME SETTINGS -------------------------"


        self.main_frame = QFrame() # Establish the main frame

        # Create the new frames and the layout
        self.main_layout = QVBoxLayout()
        self.main_layout.setContentsMargins(0, 0, 0, 0) 
        self.main_layout.setSpacing(0)
        self.up_frame = QFrame() # Upper frame
        self.bottom_frame = QFrame() # Bottom frame

        # Add the frames to the layout
        self.main_layout.addWidget(self.up_frame, 1)
        self.main_layout.addWidget(self.bottom_frame, 7)
        self.main_frame.setLayout(self.main_layout) # Establish the vertical layout

        # ------------------------- PRINCIPAL MAIN BUTTONS -------------------------"

        # Create the layout of the self.up_frame
        self.up_layout = QHBoxLayout()
        self.up_layout.setContentsMargins(30, 0, 30, 0) 
        self.up_layout.setSpacing(30)

        # Create the PushButtons
        self.button1 = CustomTriangleButton("Lig. Prep  ")
        self.button2 = CustomTriangleButton("Rec. Prep  ")
        self.button3 = CustomTriangleButton("Docking  ")
        self.button4 = CustomTriangleButton("Analysis  ")

        # Add buttons to the layout
        self.up_layout.addWidget(self.button1)
        self.up_layout.addWidget(self.button2)
        self.up_layout.addWidget(self.button3)
        self.up_layout.addWidget(self.button4)
        
        # Add the layout to the frame
        self.up_frame.setLayout(self.up_layout)

        # ------------------------- ESTABLISH BUTTON ICONS -------------------------"

        # Substrate
        substrate = QPixmap("imgs/Substrate.svg")
        substrate_icon = QIcon(substrate)
        self.button1.setIcon(substrate_icon)
        self.button1.setIconSize(substrate.size())
        self.button1.setLayoutDirection(Qt.RightToLeft)

        # Protein
        protein = QPixmap("imgs/Protein.svg")
        protein_icon = QIcon(protein)
        self.button2.setIcon(protein_icon)
        self.button2.setIconSize(protein.size())
        self.button2.setLayoutDirection(Qt.RightToLeft)

        # Docking
        docking = QPixmap("imgs/Docking.svg")
        docking_icon = QIcon(docking)
        self.button3.setIcon(docking_icon)
        self.button3.setIconSize(docking.size())
        self.button3.setLayoutDirection(Qt.RightToLeft)

        # Analysis
        analysis = QPixmap("imgs/Analysis.svg")
        analysis_icon = QIcon(analysis)
        self.button4.setIcon(analysis_icon)
        self.button4.setIconSize(analysis.size())
        self.button4.setLayoutDirection(Qt.RightToLeft)

        # ------------------------- CREATE THE MAIN STACKED WIDGET -------------------------"
        

        # Create the stacked widget
        self.main_stack = QStackedWidget()
        self.page1 = QWidget()
        self.page2 = QWidget()
        self.page3 = QWidget()
        self.page4 = QWidget()
        self.main_stack.addWidget(self.page1)
        self.main_stack.addWidget(self.page2)
        self.main_stack.addWidget(self.page3)
        self.main_stack.addWidget(self.page4)
        self.setup_page1()
        self.setup_page2()
        self.setup_page3()
        self.setup_page4()
        self.main_stack.setCurrentIndex(0)

        # Create the vertical layout of the bottom
        self.bottom_layout = QHBoxLayout()
        self.bottom_layout.setContentsMargins(0, 0, 0, 0) 
        self.bottom_layout.setSpacing(0)
        self.bottom_layout.addWidget(self.main_stack)
        

        # Establish the layout
        self.bottom_frame.setLayout(self.bottom_layout)

        # ------------------------- ESTABLISH THE STYLES -------------------------"
        self.up_frame.setStyleSheet("""
         QFrame {
                border-top: none;
                border-left: none;
                border-right: none;
                border-bottom: 2.5px solid #6d6d6d;
                background-color: #333333;
            }
        """)

        self.main_frame.setStyleSheet("""
         QFrame {
                background-color: white;
            }
        """)
        
        
        # ------------------------- ESTABLISH THE DROPDOWNWINDOW INSTANCES -------------------------"
        

        # Drop declarations
        self.dropdown_1 = DropdownWindow1(self)
        self.dropdown_1.button1Clicked.connect(self.page_11_signal)
        self.dropdown_1.button2Clicked.connect(self.page_12_signal)

        self.dropdown_2 = DropdownWindow2(self)
        self.dropdown_2.button1Clicked.connect(self.page_21_signal)


        self.dropdown_3 = DropdownWindow3(self)
        self.dropdown_3.button1Clicked.connect(self.page_31_signal)
        self.dropdown_3.button2Clicked.connect(self.page_32_signal)
        self.dropdown_3.button3Clicked.connect(self.page_33_signal)


        self.dropdown_4 = DropdownWindow4(self)
        self.dropdown_4.button1Clicked.connect(self.page_41_signal)
        self.dropdown_4.button2Clicked.connect(self.page_42_signal)
        self.dropdown_4.button3Clicked.connect(self.page_43_signal)
        self.dropdown_4.button4Clicked.connect(self.page_44_signal)
        self.dropdown_4.button5Clicked.connect(self.page_45_signal)




        # ------------------------- ESTABLISH THE MAIN WIDGET -------------------------"



        self.setCentralWidget(self.main_frame)
        


        # ------------------------- ESTABLISH THE CONNECTIONS -------------------------"
        

        self.button1.clicked.connect(self.show_dropdown_1)
        self.button2.clicked.connect(self.show_dropdown_2)
        self.button3.clicked.connect(self.show_dropdown_3)
        self.button4.clicked.connect(self.show_dropdown_4)



        # ------------------------- DROPSHOW FUNCTIONS -------------------------"


    def show_dropdown_1(self):
        if not self.ws_1:
            button_rect = self.button1.rect()
            button_pos = self.button1.mapToGlobal(button_rect.bottomLeft())
            self.dropdown_1.move(button_pos)
            self.dropdown_1.show(self.button1)
    
    
    def show_dropdown_2(self):
        if not self.ws_2:
            button_rect = self.button2.rect()
            button_pos = self.button2.mapToGlobal(button_rect.bottomLeft())
            self.dropdown_2.move(button_pos)
            self.dropdown_2.show(self.button2)

    def show_dropdown_3(self):
        if not self.ws_3:
            button_rect = self.button3.rect()
            button_pos = self.button3.mapToGlobal(button_rect.bottomLeft())
            self.dropdown_3.move(button_pos)
            self.dropdown_3.show(self.button3)
    
    def show_dropdown_4(self):
        if not self.ws_4:
            button_rect = self.button4.rect()
            button_pos = self.button4.mapToGlobal(button_rect.bottomLeft())
            self.dropdown_4.move(button_pos)
            self.dropdown_4.show(self.button4)
    
    def simulate_button_click_1(self):
        self.ws_1 = True
        self.button1.click()
        self.ws_1 = False
    
    def simulate_button_click_2(self):
        self.ws_2 = True
        self.button2.click()
        self.ws_2 = False
    
    def simulate_button_click_3(self):
        self.ws_3 = True
        self.button3.click()
        self.ws_3 = False
    
    def simulate_button_click_4(self):
        self.ws_4 = True
        self.button4.click()
        self.ws_4 = False

        # ------------------------- ESTABLISH THE PAGES FUNCTIONS -------------------------"


    # *************************************** PAGE 1 ***************************************


    def setup_page1(self):

        # Establish the stacked widget
        self.page_1_stack = QStackedWidget()
        self.page_11 = QWidget()
        self.page_12 = QWidget()
        
        # Add widgets to the stacked widget
        self.page_1_stack.addWidget(self.page_11)
        self.page_1_stack.addWidget(self.page_12)
        
        # configure the pages content
        self.setup_page11()
        self.setup_page12()
        
        # Establish the current index
        self.page_1_stack.setCurrentIndex(0)
        
        # Configure the page 1 layout
        page1_layout = QVBoxLayout()
        page1_layout.addWidget(self.page_1_stack)
        self.page1.setLayout(page1_layout)



    def setup_page11(self):

        #Principal layout 
        layout = QVBoxLayout()
        layout.setContentsMargins(20, 0, 20, 0)
        layout.setSpacing(0) 

        # Description label
        des_label = CustomTitleLabel("PREPARE LIGAND WITH MEEKO")

        # Title layout
        title_layout = QHBoxLayout()
        TSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        TSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # Label of Combo Box
        label_combo = CustomLabel("File Format: ")
        label_combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
        # ComboBox
        Combo = CustomCombo()
        Combo.addItem("smi format")
        Combo.addItem("sdf format")
        Combo.addItem("txt format")
        Combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        # File Format Spacer 2
        FFS2 = QSpacerItem(1000, 10, QSizePolicy.Expanding, QSizePolicy.Expanding)
        

        # HCombolayout
        HCombolayout = QHBoxLayout()


        # File path
        label_file = CustomLabel("File Path: ")
        label_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        label_file.setFixedHeight(30)

        # File Label for path
        label_path = CustomPathLabel("")

        # Button for file path
        Button_path = QPushButton("...")
        Button_path.setFixedWidth(40)
        Button_path.setFixedHeight(40)
        Button_path.clicked.connect(lambda: self.open_file_dialog(self.get_file_filter_11(Combo), label_path))

        # Horizontal layout
        Hlayout = QHBoxLayout()

        # Convert layout
        Clayout = QHBoxLayout()
        CSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        CSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # Line Edit
        TextLE = CustomLabel("Name folder: ")
        MyLE = CustomLineEdit()
        validator = QRegularExpressionValidator(QRegularExpression("^[A-Za-z ]*$"))
        MyLE.setValidator(validator)


        # PushButton
        PButton = CustomButton("Convert")
        PButton.clicked.connect(lambda: self.conversion_to_pdbqt_11(label_path, ProgresBar, Combo, MyLE))
        PButton.setFixedWidth(300)


        # Progress Bar
        ProgresBar = QProgressBar()
        ProgresBar.setFixedHeight(50)

        # bottom Spacer
        spacer_bottom = QSpacerItem(10, 600, QSizePolicy.Expanding, QSizePolicy.Expanding)


        #  ------------- Add Items -------------

        # Title layout
        title_layout.addItem(TSL)
        title_layout.addWidget(des_label)
        title_layout.addItem(TSR)
        layout.addItem(title_layout)

        # Horizontal combo layout
        layout.addSpacing(100)
        HCombolayout.addWidget(label_combo)
        HCombolayout.addWidget(Combo)
        HCombolayout.addItem(FFS2)
        layout.addItem(HCombolayout)
        layout.addSpacing(50)
        layout.addWidget(label_file)
        layout.setSpacing(5)


        # Horizontal label path and button path layout
        Hlayout.addWidget(label_path)
        Hlayout.addWidget(Button_path)
        layout.addItem(Hlayout)
        layout.addSpacing(50)

        # Add the line edit
        layout.addWidget(TextLE)
        layout.addWidget(MyLE)
        layout.addSpacing(50)

        # Horizontal layout for the button
        Clayout.addItem(CSL)
        Clayout.addWidget(PButton)
        Clayout.addItem(CSR)
        layout.addItem(Clayout)

        # Progress bar and bottom
        layout.addSpacing(20)
        layout.addWidget(ProgresBar)
        layout.addItem(spacer_bottom)
        
        # Set layout
        self.page_11.setLayout(layout)
    
    def conversion_to_pdbqt_11(self, label, bar, combo, LineEdit):
        
        file_path = label.text()
        folder_text = LineEdit.text() if LineEdit.text() != "" else "PDBQT files"

        if os.path.isfile(file_path):

            if (combo.currentIndex() == 0):

                try:
                    Con = Create_pdbqt_by_smi_meeko.Conversions()
                    maxim = Create_pdbqt_by_smi_meeko.Maximum()
                    maxim.contar_maximo(label.text())
                    bar.setMaximum(maxim.maximo)
                    
                    for current_count, _, _ in Con.conversion(file_path, folder_text):
                        bar.setValue(current_count)

                except Exception as e:
                    QMessageBox.critical(self, "Error 1_1-2", e)

            elif (combo.currentIndex() == 1):
                try:
                    Con = Create_pdbqt_by_sdf_meeko.Conversions()
                    maxim = Create_pdbqt_by_sdf_meeko.Maximum()

                    maxim.contar_maximo(label.text())

                    bar.setMaximum(maxim.maximo)
                    
                    for current_count, _, _ in Con.conversion(file_path, folder_text):
                        bar.setValue(current_count)

                except Exception as e:
                    QMessageBox.critical(self, "Error 1_1-3", e)

            elif (combo.currentIndex() == 2):
                try:
                    Con = Create_pdbqt_by_txt_meeko.Conversions()
                    maxim = Create_pdbqt_by_txt_meeko.Maximum()

                    maxim.contar_maximo(label.text())

                    bar.setMaximum(maxim.maximo)
                    
                    for current_count, _, _ in Con.conversion(file_path, folder_text):
                        bar.setValue(current_count)

                except Exception as e:
                    QMessageBox.critical(self, "Error 1_1-4", e)
        
        else:
            QMessageBox.critical(self, "Error 1_1-1", "Please, introduce a correct file path.")


    def get_file_filter_11(self, combo):
        if combo.currentIndex() == 0:
            return "smi files (*.smi);;All files (*)"
        elif combo.currentIndex() == 1:
            return "sdf files (*.sdf);;All files (*)"
        elif combo.currentIndex() == 2:
            return "txt files (*.txt);;All files (*)"



    def setup_page12(self):

        #  ------------- Create items -------------

        # Create the layout
        layout = QVBoxLayout()

        # Description label
        des_label = CustomTitleLabel("PREPARE LIGAND WITH OPENBABEL")

        # Title layout
        title_layout = QHBoxLayout()
        TSL = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        TSR = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # Label of Combo Box
        label_combo = CustomLabel("File Format: ")
        label_combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
        # ComboBox
        Combo = CustomCombo()
        Combo.addItem("smi format")
        Combo.addItem("sdf format")
        Combo.addItem("txt format")
        Combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        FFS2 = QSpacerItem(1000, 10, QSizePolicy.Expanding, QSizePolicy.Expanding)

        # File path
        label_file = CustomLabel("File Path: ")
        label_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        label_file.setFixedHeight(30)

        # File Label for path
        label_path = CustomPathLabel("")

        # Button for file path
        Button_path = QPushButton("...")
        Button_path.setFixedWidth(40)
        Button_path.setFixedHeight(40)
        Button_path.clicked.connect(lambda: self.open_file_dialog(self.get_file_filter_12(Combo), label_path))

        # Horizontal layout
        Hlayout = QHBoxLayout()

        # Line Edit
        TextLE = CustomLabel("Name folder: ")
        MyLE = CustomLineEdit()
        validator = QRegularExpressionValidator(QRegularExpression("^[A-Za-z ]*$"))
        MyLE.setValidator(validator)

        # PushButton
        PButton = CustomButton("Convert")
        PButton.clicked.connect(lambda: self.conversion_to_pdbqt_12(label_path, ProgresBar, Combo, MyLE))
        PButton.setFixedWidth(300)

        # Progress Bar
        ProgresBar = QProgressBar()
        ProgresBar.setFixedHeight(50)

        # Convert layout
        Clayout = QHBoxLayout()
        CSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        CSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        
        # Bottom Spacer
        spacer_bottom = QSpacerItem(20, 10, QSizePolicy.Minimum, QSizePolicy.Expanding)


        #  ------------- Add Items -------------
        
        # Add title
        title_layout.addItem(TSL)
        title_layout.addWidget(des_label)
        title_layout.addItem(TSR)
        layout.addItem(title_layout)

        # Add Combo
        HCombolayout = QHBoxLayout()
        layout.addSpacing(100)
        HCombolayout.addWidget(label_combo)
        HCombolayout.addWidget(Combo)
        HCombolayout.addItem(FFS2)
        layout.addItem(HCombolayout)
        layout.addSpacing(50)

        # Horizontal label path and button path layout
        layout.addWidget(label_file)
        Hlayout.addWidget(label_path)
        Hlayout.addWidget(Button_path)
        layout.addItem(Hlayout)
        layout.addSpacing(50)

        # Add the line edit
        layout.addWidget(TextLE)
        layout.addWidget(MyLE)
        layout.addSpacing(50)
        
        # Horizontal layout for the button
        Clayout.addItem(CSL)
        Clayout.addWidget(PButton)
        Clayout.addItem(CSR)
        layout.addItem(Clayout)

        # Progress bar and bottom
        layout.addSpacing(20)
        layout.addWidget(ProgresBar)
        layout.addItem(spacer_bottom)

        
        # Set layout and final spacer
        layout.addItem(spacer_bottom)
        self.page_12.setLayout(layout)

    def conversion_to_pdbqt_12(self, label, bar, combo, LineEdit):
        
        file_path = label.text()
        folder_text = LineEdit.text() if LineEdit.text() != "" else "PDBQT files"

        if os.path.isfile(file_path):

            if (combo.currentIndex() == 0):

                try:
                    Con = Create_pdbqt_by_smi_obabel.Conversions()
                    Con.Maximum(label.text())
                    bar.setMaximum(Con.maxim)

                    for pdbqt_file in Con.conversions(folder_text, file_path):
                        bar.setValue(Con.contator)

                except Exception as e:
                    QMessageBox.critical(self, "Error 1_2-1", e)

            elif (combo.currentIndex() == 1):
                try:
                    Con = Create_pdbqt_by_sdf_obabel.Conversions()
                    Con.Maximum(label.text())
                    bar.setMaximum(Con.maxim)
                    
                    
                    for pdbqt_file in Con.conversions(folder_text, file_path):
                        bar.setValue(Con.contator)
                    

                except Exception as e:
                    QMessageBox.critical(self, "Error 1_2-2", e)

            elif (combo.currentIndex() == 2):
                
                try:

                    Con = Create_pdbqt_by_txt_obabel.Conversions()
                    Con.Maximum(label.text())
                    bar.setMaximum(Con.maxim)

                    for pdbqt_file in Con.conversions(folder_text, file_path):
                        bar.setValue(Con.contator)
                
                except Exception as e:
                    QMessageBox.critical(self, "Error 1_2-3", e)

        
        else:
            QMessageBox.critical(self, "Error 1_2-1", "Please, introduce a correct file path.")

    def get_file_filter_12(self, combo):
        if combo.currentIndex() == 0:
            return "smi files (*.smi);;All files (*)"
        elif combo.currentIndex() == 1:
            return "sdf files (*.sdf);;All files (*)"
        elif combo.currentIndex() == 2:
            return "txt files (*.txt);;All files (*)"


    # *************************************** PAGE 2 ***************************************

    def setup_page2(self):
        
        # Establish the stacked widget
        self.page_2_stack = QStackedWidget()
        self.page_21 = QWidget()
        self.page_22 = QWidget()
        
        # Add widgets to the stacked widget
        self.page_2_stack.addWidget(self.page_21)
        self.page_2_stack.addWidget(self.page_22)
        
        # configure the pages content
        self.setup_page21()
        
        # Establish the current index
        self.page_2_stack.setCurrentIndex(0)
        
        # Configure the page 2 layout
        page2_layout = QVBoxLayout()
        page2_layout.addWidget(self.page_2_stack)
        self.page2.setLayout(page2_layout)


    def setup_page21(self):

        #  ------------- Create items -------------

        # Create the layout
        layout = QVBoxLayout()

        # Description label
        des_label = CustomTitleLabel("PREPARE RECEPTOR WITH OPENBABEL")

        # Title layout
        title_layout = QHBoxLayout()
        TSL = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        TSR = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        
        # File path
        label_file = CustomLabel("Folder Path: ")
        label_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        label_file.setFixedHeight(30)

        # File Label for path
        label_path = CustomPathLabel("")

        # Button for file path
        Button_path = QPushButton("...")
        Button_path.setFixedWidth(40)
        Button_path.setFixedHeight(40)
        Button_path.clicked.connect(lambda: self.open_folder_dialog(label_path))

        # Horizontal layout
        Hlayout = QHBoxLayout()

        # Line Edit
        TextLE = CustomLabel("Name folder: ")
        MyLE = CustomLineEdit()
        validator = QRegularExpressionValidator(QRegularExpression("^[A-Za-z ]*$"))
        MyLE.setValidator(validator)

        # Checkbox and his text
        Checklayout = QHBoxLayout()
        TextCheck = CustomLabel("Save ligands: ")
        MyCheck = CustomCheckBox()
        Checkspacer = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # PushButton
        PButton = CustomButton("Convert")
        PButton.clicked.connect(lambda: self.conversion_to_pdbqt_21(label_path, ProgresBar, MyLE, MyCheck))
        PButton.setFixedWidth(300)

        # Progress Bar
        ProgresBar = QProgressBar()
        ProgresBar.setFixedHeight(50)

        # Convert layout
        Clayout = QHBoxLayout()
        CSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        CSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        
        # Bottom Spacer
        spacer_bottom = QSpacerItem(20, 10, QSizePolicy.Minimum, QSizePolicy.Expanding)

        # Create layout
        layout = QVBoxLayout()

        #  ------------- Add Items -------------
        
        # Add title
        title_layout.addItem(TSL)
        title_layout.addWidget(des_label)
        title_layout.addItem(TSR)
        layout.addItem(title_layout)
        layout.addSpacing(50)
        layout.addWidget(label_file)

        # Horizontal label path and button path layout
        Hlayout.addWidget(label_path)
        Hlayout.addWidget(Button_path)
        layout.addItem(Hlayout)
        layout.addSpacing(50)

        # Add the line edit
        layout.addWidget(TextLE)
        layout.addWidget(MyLE)
        layout.addSpacing(50)

        # Add the checkbox
        Checklayout.addWidget(TextCheck)
        Checklayout.addWidget(MyCheck)
        Checklayout.addItem(Checkspacer)
        layout.addItem(Checklayout)
        layout.addSpacing(50)
        
        # Horizontal layout for the button
        Clayout.addItem(CSL)
        Clayout.addWidget(PButton)
        Clayout.addItem(CSR)
        layout.addItem(Clayout)

        # Progress bar and bottom
        layout.addSpacing(20)
        layout.addWidget(ProgresBar)
        layout.addItem(spacer_bottom)

        # Set layout and final spacer
        layout.addItem(spacer_bottom)
        self.page_12.setLayout(layout)

        # Establish layout
        self.page_21.setLayout(layout)

    def conversion_to_pdbqt_21(self, label, bar, LineEdit, CheckBox):
        
        file_path = label.text()
        folder_text = LineEdit.text() if LineEdit.text() != "" else "PDBQT files"
        status = CheckBox.isChecked()

        if os.path.isdir(file_path):

            try:
                Con = Create_protein_pdbqt_by_obabel.Conversions()
                Con.Maximum(file_path)
                bar.setMaximum(Con.maxim)
                list_chains = []
                
                Con.get_proteins(file_path)
                list_proteins = Con.proteins

                for pdb_file in Con.Chains(file_path):
                    list_chains.append(Con.chains)
                
                
                generator = Con.conversions(file_path, folder_text, status)
                for i, chain in enumerate(list_chains):
                    if "B" in chain:
                        dialog = CheckableOptionsDialog(chain, list_proteins[i], self)

                        if dialog.exec():
                            selected_options = dialog.get_checked_options()
                            Con.selected_chains = selected_options
                    else:
                        Con.selected_chains = chain
                    
                    next(generator)
                    bar.setValue(Con.contator)
                

            except Exception as e:
                QMessageBox.critical(self, "Error 2_1-2", e)

        
        else:
            QMessageBox.critical(self, "Error 2_1-1", "Please, introduce a correct file path.")


    # *************************************** PAGE 3 ***************************************

    def setup_page3(self):

        #  ------------- Create items -------------

        # Establish the stacked widget
        self.page_3_stack = QStackedWidget()
        self.page_31 = QWidget()
        self.page_32 = QWidget()
        self.page_33 = QWidget()
        
        # Add widgets to the stacked widget
        self.page_3_stack.addWidget(self.page_31)
        self.page_3_stack.addWidget(self.page_32)
        self.page_3_stack.addWidget(self.page_33)
        
        # configure the pages content
        self.setup_page31()
        self.setup_page32()
        self.setup_page33()
        
        # Establish the current index
        self.page_3_stack.setCurrentIndex(0)
        
        # Configure the page 2 layout
        page3_layout = QVBoxLayout()
        page3_layout.addWidget(self.page_3_stack)
        self.page3.setLayout(page3_layout)

    def setup_page31(self):

        
        # Create the layout
        layout = QVBoxLayout()

        # Description label
        des_label = CustomTitleLabel("MOLECULAR DOCKING WITH SMINA")

        # Title layout
        title_layout = QHBoxLayout()
        TSL = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        TSR = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        
        # Select protein
        # File path
        label_file = CustomLabel("Protein path: ")
        label_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        label_file.setFixedHeight(30)

        # File Label for path
        label_path = CustomPathLabel("")

        # Button for file path
        Button_path = QPushButton("...")
        Button_path.setFixedWidth(40)
        Button_path.setFixedHeight(40)
        Button_path.clicked.connect(lambda: self.open_file_dialog("pdbqt files (*.pdbqt);;All files (*)", label_path))

        # Horizontal layout
        Hlayout = QHBoxLayout()

        # File path
        labelligand_file = CustomLabel("Ligand folder: ")
        labelligand_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        labelligand_file.setFixedHeight(30)

        # File label for ligands path
        ligands_path = CustomPathLabel("")

        # Button for ligands path
        Button_lpath = QPushButton("...")
        Button_lpath.setFixedWidth(40)
        Button_lpath.setFixedHeight(40)
        Button_lpath.clicked.connect(lambda: self.open_folder_dialog(ligands_path))

        # Horizontal ligand layout
        hligandlayout = QHBoxLayout()

        # Name folder
        # Line Edit
        TextLE = CustomLabel("Name folder: ")
        MyLE = CustomLineEdit()
        validator = QRegularExpressionValidator(QRegularExpression("^[A-Za-z ]*$"))
        MyLE.setValidator(validator)

        # Horizontal size layout
        slayout = QHBoxLayout()

        # Validator for numbers
        size_validator = QRegularExpressionValidator(QRegularExpression("^-?[0-9]*\\.?[0-9]+$"))

        # Labels for size and Coordinates
        XLab = CustomLabel("X: ")
        YLab = CustomLabel("Y: ")
        ZLab = CustomLabel("Z: ")
        X_SLab = CustomLabel("sX: ")
        Y_SLab = CustomLabel("sY: ")
        Z_SLab = CustomLabel("sZ: ")

        # Line edits for size and Coordinates
        XLE = CustomLineEdit()
        XLE.setValidator(size_validator)
        YLE = CustomLineEdit()
        YLE.setValidator(size_validator)
        ZLE = CustomLineEdit()
        ZLE.setValidator(size_validator)
        X_SLE = CustomLineEdit()
        X_SLE.setValidator(size_validator)
        Y_SLE = CustomLineEdit()
        Y_SLE.setValidator(size_validator)
        Z_SLE = CustomLineEdit()
        Z_SLE.setValidator(size_validator)

        # Label of Combo Box
        label_combo = CustomLabel("Scoring: ")
        label_combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
        # ComboBox
        Combo = CustomCombo()
        Combo.addItem("Vina")
        Combo.addItem("Vinardo")
        Combo.addItem("Dkoes")
        Combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        FFS2 = QSpacerItem(1000, 10, QSizePolicy.Expanding, QSizePolicy.Expanding)

        # PushButton
        PButton = CustomButton("Do docking")
        PButton.clicked.connect(lambda: self.conversion_to_pdbqt_31(label_path, ligands_path, MyLE, XLE, YLE, ZLE, X_SLE, Y_SLE, Z_SLE, ProgresBar, Combo))
        PButton.setFixedWidth(300)

        # Progress Bar
        ProgresBar = QProgressBar()
        ProgresBar.setFixedHeight(50)

        # Convert layout
        Clayout = QHBoxLayout()
        CSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        CSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # Bottom Spacer
        spacer_bottom = QSpacerItem(20, 10, QSizePolicy.Minimum, QSizePolicy.Expanding)

        #  ------------------------- Add Items -------------------------

        # Add the title
        title_layout.addItem(TSL)
        title_layout.addWidget(des_label)
        title_layout.addItem(TSR)
        layout.addItem(title_layout)
        layout.addSpacing(10)

        # Add protein label path
        layout.addWidget(label_file)
        Hlayout.addWidget(label_path)
        Hlayout.addWidget(Button_path)
        layout.addItem(Hlayout)
        layout.addSpacing(10)

        # Add ligand label path
        layout.addWidget(labelligand_file)
        hligandlayout.addWidget(ligands_path)
        hligandlayout.addWidget(Button_lpath)
        layout.addItem(hligandlayout)
        layout.addSpacing(10)

        # Add the line edit
        layout.addWidget(TextLE)
        layout.addWidget(MyLE)
        layout.addSpacing(30)

        # Add the sizes
        slayout.addWidget(XLab)
        slayout.addWidget(XLE)
        slayout.addSpacing(30)
        slayout.addWidget(YLab)
        slayout.addWidget(YLE)
        slayout.addSpacing(30)
        slayout.addWidget(ZLab)
        slayout.addWidget(ZLE)
        slayout.addSpacing(30)
        slayout.addWidget(X_SLab)
        slayout.addWidget(X_SLE)
        slayout.addSpacing(30)
        slayout.addWidget(Y_SLab)
        slayout.addWidget(Y_SLE)
        slayout.addSpacing(30)
        slayout.addWidget(Z_SLab)
        slayout.addWidget(Z_SLE)
        layout.addItem(slayout)
        layout.addSpacing(20)

        # Add Combo
        HCombolayout = QHBoxLayout()
        HCombolayout.addWidget(label_combo)
        HCombolayout.addWidget(Combo)
        HCombolayout.addItem(FFS2)
        layout.addItem(HCombolayout)


        # Horizontal layout for the button
        Clayout.addItem(CSL)
        Clayout.addWidget(PButton)
        Clayout.addItem(CSR)
        layout.addItem(Clayout)

        # Progress bar and bottom
        layout.addSpacing(20)
        layout.addWidget(ProgresBar)
        layout.addItem(spacer_bottom)

        # Add bottom spacer
        layout.addItem(spacer_bottom)

        # Set layout
        self.page_31.setLayout(layout) 

    def conversion_to_pdbqt_31(self, protein_label, ligand_label, folder_LineEdit, X_LE, Y_LE, Z_LE, XS_LE, YS_LE, ZS_LE, bar, combo):
        
        # Recovery the text
        protein_path = protein_label.text()
        ligand_path = ligand_label.text()
        folder_text = folder_LineEdit.text() if folder_LineEdit.text() != "" else "Docking files"
        X = X_LE.text()
        Y = Y_LE.text()
        Z = Z_LE.text()
        XS = XS_LE.text()
        YS = YS_LE.text()
        ZS = ZS_LE.text()
        combo_index = combo.currentIndex()

        if os.path.isfile(protein_path) and os.path.isdir(ligand_path) and X and Y and Z and XS and YS and ZS:
            
            try:
                Con = Docking_with_smina.Conversions()
                Con.Maximum(ligand_path)
                bar.setMaximum(Con.maxim)

                for ligand in Con.conversions(protein_path, ligand_path, folder_text, X, Y, Z, XS, YS, ZS, combo_index):
                    bar.setValue(Con.contator)

            except Exception as e:
                QMessageBox.critical(self, "Error 3_1-2", e)
            
        else:
            QMessageBox.critical(self, "Error 3_1-1", "Make sure to fill in all fields.")


    def setup_page32(self):

        # Create the layout
        layout = QVBoxLayout()

        # Description label
        des_label = CustomTitleLabel("REDOCKING WITH SMINA")

        # Title layout
        title_layout = QHBoxLayout()
        TSL = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        TSR = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        
        # Select protein
        # File path
        label_file = CustomLabel("Protein path: ")
        label_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        label_file.setFixedHeight(30)

        # File Label for path
        label_path = CustomPathLabel("")

        # Button for file path
        Button_path = QPushButton("...")
        Button_path.setFixedWidth(40)
        Button_path.setFixedHeight(40)
        Button_path.clicked.connect(lambda: self.open_file_dialog("pdbqt files (*.pdbqt);;All files (*)", label_path))

        # Horizontal layout
        Hlayout = QHBoxLayout()

        # File path
        labelligand_file = CustomLabel("Ligand path: ")
        labelligand_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        labelligand_file.setFixedHeight(30)

        # File label for ligands path
        ligands_path = CustomPathLabel("")

        # Button for ligands path
        Button_lpath = QPushButton("...")
        Button_lpath.setFixedWidth(40)
        Button_lpath.setFixedHeight(40)
        Button_lpath.clicked.connect(lambda: self.open_file_dialog("pdbqt files (*.pdbqt);;All files (*)", ligands_path))

        # Horizontal ligand layout
        hligandlayout = QHBoxLayout()

        # Name folder
        # Line Edit
        TextLE = CustomLabel("Name folder: ")
        MyLE = CustomLineEdit()
        validator = QRegularExpressionValidator(QRegularExpression("^[A-Za-z ]*$"))
        MyLE.setValidator(validator)

        # Horizontal size layout
        slayout = QHBoxLayout()

        # Validator for numbers
        size_validator = QRegularExpressionValidator(QRegularExpression("^-?[0-9]*\\.?[0-9]+$"))


        # Label of Combo Box
        label_combo = CustomLabel("Scoring: ")
        label_combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
        # ComboBox
        Combo = CustomCombo()
        Combo.addItem("Vina")
        Combo.addItem("Vinardo")
        Combo.addItem("Dkoes")
        Combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        FFS2 = QSpacerItem(1000, 10, QSizePolicy.Expanding, QSizePolicy.Expanding)

        # PushButton
        PButton = CustomButton("Do docking")
        PButton.clicked.connect(lambda: self.conversion_to_pdbqt_32(label_path, ligands_path, MyLE, ProgresBar, Combo))
        PButton.setFixedWidth(300)

        # Progress Bar
        ProgresBar = QProgressBar()
        ProgresBar.setFixedHeight(50)

        # Convert layout
        Clayout = QHBoxLayout()
        CSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        CSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # Bottom Spacer
        spacer_bottom = QSpacerItem(20, 10, QSizePolicy.Minimum, QSizePolicy.Expanding)

        #  ------------------------- Add Items -------------------------

        # Add the title
        title_layout.addItem(TSL)
        title_layout.addWidget(des_label)
        title_layout.addItem(TSR)
        layout.addItem(title_layout)
        layout.addSpacing(20)

        # Add protein label path
        layout.addWidget(label_file)
        Hlayout.addWidget(label_path)
        Hlayout.addWidget(Button_path)
        layout.addItem(Hlayout)
        layout.addSpacing(20)

        # Add ligand label path
        layout.addWidget(labelligand_file)
        hligandlayout.addWidget(ligands_path)
        hligandlayout.addWidget(Button_lpath)
        layout.addItem(hligandlayout)
        layout.addSpacing(20)

        # Add the line edit
        layout.addWidget(TextLE)
        layout.addWidget(MyLE)
        layout.addSpacing(40)

        # Add Combo
        HCombolayout = QHBoxLayout()
        HCombolayout.addWidget(label_combo)
        HCombolayout.addWidget(Combo)
        HCombolayout.addItem(FFS2)
        layout.addItem(HCombolayout)


        # Horizontal layout for the button
        Clayout.addItem(CSL)
        Clayout.addWidget(PButton)
        Clayout.addItem(CSR)
        layout.addItem(Clayout)

        # Progress bar and bottom
        layout.addSpacing(20)
        layout.addWidget(ProgresBar)
        layout.addItem(spacer_bottom)

        # Add bottom spacer
        layout.addItem(spacer_bottom)

        # Set layout
        self.page_32.setLayout(layout)

    def conversion_to_pdbqt_32(self, protein_label, ligand_label, folder_LineEdit, bar, combo):
        
        # Recovery the text
        protein_path = protein_label.text()
        ligand_path = ligand_label.text()
        folder_text = folder_LineEdit.text() if folder_LineEdit.text() != "" else "Redocking files"
        combo_index = combo.currentIndex()

        if os.path.isfile(protein_path) and os.path.isfile(ligand_path):
            
            try:
                Con = Redocking_with_smina.Conversions()
                Con.Maximum(ligand_path)
                bar.setMaximum(Con.maxim)

                Con.conversions(protein_path, ligand_path, folder_text, combo_index)
                bar.setValue(Con.contator)

            except Exception as e:
                QMessageBox.critical(self, "Error 3_1-2", e)
            
        else:
            QMessageBox.critical(self, "Error 3_1-1", "Make sure to fill in all fields.")

    def setup_page33(self):

        layout = QVBoxLayout()
        label = CustomLabel("Este es un QLabel en la Página 33")
        CB = CustomCheckBox()
        button = QPushButton("Este es un QPushButton en la Página 33") 
        layout.addWidget(label)
        layout.addWidget(CB)
        layout.addWidget(button)
        self.page_33.setLayout(layout) 



    # *************************************** PAGE 4 ***************************************

    def setup_page4(self):

        # Establish the stacked widget
        self.page_4_stack = QStackedWidget()
        self.page_41 = QWidget()
        self.page_42 = QWidget()
        self.page_43 = QWidget()
        self.page_44 = QWidget()
        self.page_45 = QWidget()
        self.page_46 = QWidget()
        self.page_47 = QWidget()
        self.page_48 = QWidget()
        
        # Add widgets to the stacked widget
        self.page_4_stack.addWidget(self.page_41)
        self.page_4_stack.addWidget(self.page_42)
        self.page_4_stack.addWidget(self.page_43)
        self.page_4_stack.addWidget(self.page_44)
        self.page_4_stack.addWidget(self.page_45)
        self.page_4_stack.addWidget(self.page_46)
        self.page_4_stack.addWidget(self.page_47)
        self.page_4_stack.addWidget(self.page_48)
        
        # configure the pages content
        self.setup_page41()
        self.setup_page42()
        self.setup_page43()
        self.setup_page44()
        self.setup_page45()
        self.setup_page46()
        self.setup_page47()
        self.setup_page48()
        
        # Establish the current index
        self.page_4_stack.setCurrentIndex(0)
        
        # Configure the page 4 layout
        page4_layout = QVBoxLayout()
        page4_layout.addWidget(self.page_4_stack)
        self.page4.setLayout(page4_layout)


    def setup_page41(self):

        #Principal layout 
        layout = QVBoxLayout()
        layout.setContentsMargins(20, 0, 20, 0)
        layout.setSpacing(0) 

        # Description label
        des_label = CustomTitleLabel("SPLIT PDBQT FILES")

        # Title layout
        title_layout = QHBoxLayout()
        TSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        TSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)


        # File Format Spacer 2
        FFS2 = QSpacerItem(1000, 10, QSizePolicy.Expanding, QSizePolicy.Expanding)
        

        # HCombolayout
        HCombolayout = QHBoxLayout()


        # File path
        label_file = CustomLabel("File Path: ")
        label_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        label_file.setFixedHeight(30)

        # File Label for path
        label_path = CustomPathLabel("")

        # Button for file path
        Button_path = QPushButton("...")
        Button_path.setFixedWidth(40)
        Button_path.setFixedHeight(40)
        Button_path.clicked.connect(lambda: self.open_folder_dialog(label_path))

        # Horizontal layout
        Hlayout = QHBoxLayout()

        # Convert layout
        Clayout = QHBoxLayout()
        CSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        CSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # Line Edit
        TextLE = CustomLabel("Name folder: ")
        MyLE = CustomLineEdit()
        validator = QRegularExpressionValidator(QRegularExpression("^[A-Za-z ]*$"))
        MyLE.setValidator(validator)


        # PushButton
        PButton = CustomButton("Convert")
        PButton.clicked.connect(lambda: self.conversion_to_pdbqt_41(label_path, MyLE , ProgresBar))
        PButton.setFixedWidth(300)


        # Progress Bar
        ProgresBar = QProgressBar()
        ProgresBar.setFixedHeight(50)

        # bottom Spacer
        spacer_bottom = QSpacerItem(10, 600, QSizePolicy.Expanding, QSizePolicy.Expanding)

        #  ------------- Add Items -------------

        # Title layout
        title_layout.addItem(TSL)
        title_layout.addWidget(des_label)
        title_layout.addItem(TSR)
        layout.addItem(title_layout)

        # Horizontal combo layout
        layout.addSpacing(100)
        layout.addItem(HCombolayout)
        layout.addSpacing(50)
        layout.addWidget(label_file)
        layout.setSpacing(5)


        # Horizontal label path and button path layout
        Hlayout.addWidget(label_path)
        Hlayout.addWidget(Button_path)
        layout.addItem(Hlayout)
        layout.addSpacing(50)

        # Add the line edit
        layout.addWidget(TextLE)
        layout.addWidget(MyLE)
        layout.addSpacing(50)

        # Horizontal layout for the button
        Clayout.addItem(CSL)
        Clayout.addWidget(PButton)
        Clayout.addItem(CSR)
        layout.addItem(Clayout)

        # Progress bar and bottom
        layout.addSpacing(20)
        layout.addWidget(ProgresBar)
        layout.addItem(spacer_bottom)
        
        # Set layout
        self.page_41.setLayout(layout)


    def conversion_to_pdbqt_41(self, file_folder, folder_LineEdit, bar):
        
        # Recovery the text
        ligands = file_folder.text()
        folder_text = folder_LineEdit.text() if folder_LineEdit.text() != "" else "Split files"

        if os.path.isdir(ligands):
            
            try:
                Con = Convert_pdbqt_to_mol2.Conversions()
                Con.Maximum(ligands)
                bar.setMaximum(Con.maxim)

                for ligand in Con.conversions(ligands, folder_text):
                    bar.setValue(Con.contator)


            except Exception as e:
                QMessageBox.critical(self, "Error 4_1-2", e)
            
        else:
            QMessageBox.critical(self, "Error 4_1-1", "Please enter a correct path")


    def setup_page42(self):

        # Create the layout
        layout = QVBoxLayout()

        # Description label
        des_label = CustomTitleLabel("RMSD WITH DockRMSD")

        # Title layout
        title_layout = QHBoxLayout()
        TSL = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        TSR = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        
        # Select protein
        # File path
        label_file = CustomLabel("File 1 Path: ")
        label_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        label_file.setFixedHeight(30)

        # File Label for path
        label_path = CustomPathLabel("")

        # Button for file path
        Button_path = QPushButton("...")
        Button_path.setFixedWidth(40)
        Button_path.setFixedHeight(40)
        Button_path.clicked.connect(lambda: self.open_file_dialog("pdbqt files (*.pdbqt);;All files (*)", label_path))

        # Horizontal layout
        Hlayout = QHBoxLayout()

        # File path
        labelligand_file = CustomLabel("File 2 path: ")
        labelligand_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        labelligand_file.setFixedHeight(30)

        # File label for ligands path
        ligands_path = CustomPathLabel("")

        # Button for ligands path
        Button_lpath = QPushButton("...")
        Button_lpath.setFixedWidth(40)
        Button_lpath.setFixedHeight(40)
        Button_lpath.clicked.connect(lambda: self.open_file_dialog("pdbqt files (*.pdbqt);;All files (*)", ligands_path))

        # Horizontal ligand layout
        hligandlayout = QHBoxLayout()

        # Name folder
        # Line Edit
        TextLE = CustomLabel("Name folder: ")
        MyLE = CustomLineEdit()
        validator = QRegularExpressionValidator(QRegularExpression("^[A-Za-z ]*$"))
        MyLE.setValidator(validator)

        # Horizontal size layout
        slayout = QHBoxLayout()

        # Validator for numbers
        size_validator = QRegularExpressionValidator(QRegularExpression("^-?[0-9]*\\.?[0-9]+$"))


        # Label of Combo Box
        label_combo = CustomLabel("Scoring: ")
        label_combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        

        # PushButton
        PButton = CustomButton("Convert")
        PButton.clicked.connect(lambda: self.conversion_to_pdbqt_42(label_path, ligands_path, MyLE, ProgresBar))
        PButton.setFixedWidth(300)

        # Progress Bar
        ProgresBar = QProgressBar()
        ProgresBar.setFixedHeight(50)

        # Convert layout
        Clayout = QHBoxLayout()
        CSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        CSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # Bottom Spacer
        spacer_bottom = QSpacerItem(20, 10, QSizePolicy.Minimum, QSizePolicy.Expanding)

        #  ------------------------- Add Items -------------------------

        # Add the title
        title_layout.addItem(TSL)
        title_layout.addWidget(des_label)
        title_layout.addItem(TSR)
        layout.addItem(title_layout)
        layout.addSpacing(30)

        # Add protein label path
        layout.addWidget(label_file)
        Hlayout.addWidget(label_path)
        Hlayout.addWidget(Button_path)
        layout.addItem(Hlayout)
        layout.addSpacing(30)

        # Add ligand label path
        layout.addWidget(labelligand_file)
        hligandlayout.addWidget(ligands_path)
        hligandlayout.addWidget(Button_lpath)
        layout.addItem(hligandlayout)
        layout.addSpacing(30)

        # Add the line edit
        layout.addWidget(TextLE)
        layout.addWidget(MyLE)
        layout.addSpacing(50)

        # Horizontal layout for the button
        Clayout.addItem(CSL)
        Clayout.addWidget(PButton)
        Clayout.addItem(CSR)
        layout.addItem(Clayout)

        # Progress bar and bottom
        layout.addSpacing(20)
        layout.addWidget(ProgresBar)
        layout.addItem(spacer_bottom)

        # Add bottom spacer
        layout.addItem(spacer_bottom)

        # Set layout
        self.page_42.setLayout(layout)
    
    def conversion_to_pdbqt_42(self, file_1, file_2, folder_LineEdit, bar):
        
        # Recovery the text
        file_1 = file_1.text()
        file_2 = file_2.text()
        folder_text = folder_LineEdit.text() if folder_LineEdit.text() != "" else "RMSD file"

        if os.path.isfile(file_1) and os.path.isfile(file_2):
            
            try:
                Con = Calculate_RMSD.Conversions()
                Con.Maximum()
                bar.setMaximum(Con.maxim)

                Con.conversions(file_1, file_2, folder_text)
                bar.setValue(Con.contator)


            except Exception as e:
                QMessageBox.critical(self, "Error 4_1-2", e)
            
        else:
            QMessageBox.critical(self, "Error 4_2-2", "Please enter a correct path")
    
    def setup_page43(self):

        #  ------------- Create items -------------

        # Create the layout
        layout = QVBoxLayout()

        # Description label
        des_label = CustomTitleLabel("BEST COMPOUNDS")

        # Title layout
        title_layout = QHBoxLayout()
        TSL = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        TSR = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # Label of Combo Box
        label_combo = CustomLabel("Filter method: ")
        label_combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
        # ComboBox
        Combo = CustomCombo()
        Combo.addItem("Best compounds")
        Combo.addItem("Cutting Interval")
        Combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        Combo.currentIndexChanged.connect(lambda: self.Index_changed_43(NumberLA, NumberLE, Combo))
        FFS2 = QSpacerItem(1000, 10, QSizePolicy.Expanding, QSizePolicy.Expanding)

        # Line Edit for numbers
        NumberLA = CustomLabel("Amount of compounds: ")
        NumberLE = CustomLineEdit()
        NumberLE.setFixedWidth(200)
        validator = QIntValidator(0, 999999)
        NLES = QSpacerItem(100, 10, QSizePolicy.Expanding, QSizePolicy.Expanding)
        NumberLE.setValidator(validator)

        # File path
        label_file = CustomLabel("File folder: ")
        label_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        label_file.setFixedHeight(30)

        # File Label for path
        label_path = CustomPathLabel("")

        # Button for file path
        Button_path = QPushButton("...")
        Button_path.setFixedWidth(40)
        Button_path.setFixedHeight(40)
        Button_path.clicked.connect(lambda: self.open_folder_dialog(label_path))

        # Horizontal layout
        Hlayout = QHBoxLayout()

        # Line Edit
        TextLE = CustomLabel("Name folder: ")
        MyLE = CustomLineEdit()
        validator = QRegularExpressionValidator(QRegularExpression("^[A-Za-z ]*$"))
        MyLE.setValidator(validator)

        # PushButton
        PButton = CustomButton("Convert")
        PButton.clicked.connect(lambda: self.conversion_to_pdbqt_43(NumberLE, label_path, MyLE, Combo, ProgresBar))
        PButton.setFixedWidth(300)

        # Progress Bar
        ProgresBar = QProgressBar()
        ProgresBar.setFixedHeight(50)

        # Convert layout
        Clayout = QHBoxLayout()
        CSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        CSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        
        # Bottom Spacer
        spacer_bottom = QSpacerItem(20, 10, QSizePolicy.Minimum, QSizePolicy.Expanding)


        #  ------------- Add Items -------------
        
        # Add title
        title_layout.addItem(TSL)
        title_layout.addWidget(des_label)
        title_layout.addItem(TSR)
        layout.addItem(title_layout)

        # Add Combo
        HCombolayout = QHBoxLayout()
        layout.addSpacing(20)
        HCombolayout.addWidget(label_combo)
        HCombolayout.addWidget(Combo)
        HCombolayout.addItem(FFS2)
        layout.addItem(HCombolayout)
        layout.addSpacing(30)

        # Add number LE and LA
        HComboNlayout = QHBoxLayout()
        HComboNlayout.addWidget(NumberLA)
        HComboNlayout.addWidget(NumberLE)
        HComboNlayout.addItem(NLES)
        layout.addItem(HComboNlayout)
        layout.addSpacing(30)


        # Horizontal label path and button path layout
        layout.addWidget(label_file)
        Hlayout.addWidget(label_path)
        Hlayout.addWidget(Button_path)
        layout.addItem(Hlayout)
        layout.addSpacing(30)


        # Add the line edit
        layout.addWidget(TextLE)
        layout.addWidget(MyLE)
        layout.addSpacing(30)
        
        # Horizontal layout for the button
        Clayout.addItem(CSL)
        Clayout.addWidget(PButton)
        Clayout.addItem(CSR)
        layout.addItem(Clayout)

        # Progress bar and bottom
        layout.addSpacing(20)
        layout.addWidget(ProgresBar)
        layout.addItem(spacer_bottom)

        
        # Set layout and final spacer
        layout.addItem(spacer_bottom)
        self.page_43.setLayout(layout)
    
    def conversion_to_pdbqt_43(self, NLE, logs_folder, des_folder, Combo, bar):
        
        # Recovery the text
        NLE = NLE.text()
        logs = logs_folder.text()
        folder_text = des_folder.text() if des_folder.text() != "" else "Screening folder"
        index = Combo.currentIndex()

        if os.path.isdir(logs) and NLE != "":
            
            try:
                Con = Best_screening.Conversions()
                Con.Maximum(logs)
                bar.setMaximum(Con.maxim)

                for log in Con.conversions(NLE, logs, folder_text, index):
                    bar.setValue(Con.contator)

            except Exception as e:
                QMessageBox.critical(self, "Error 4_3-2", e)
            
        else:
            QMessageBox.critical(self, "Error 4_3-1", "Please fill in all fields")
    
    def Index_changed_43(self, LA, LE, combo):
        
        if combo.currentIndex() == 0:
            LE.setText("0")
            LA.setText("Amount of compunds: ")
            Validator = QIntValidator(0, 999999)
            LE.setValidator(Validator)
        elif combo.currentIndex() == 1:
            LE.setText("0")
            LA.setText("Cutting interval: ")
            Validator = QDoubleValidator()
            LE.setValidator(Validator)
    
    def setup_page44(self):

        #  ------------- Create items -------------

        # Create the layout
        layout = QVBoxLayout()

        # Description label
        des_label = CustomTitleLabel("PREDICT BBB")

        # Title layout
        title_layout = QHBoxLayout()
        TSL = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)
        TSR = QSpacerItem(10000, 10, QSizePolicy.Maximum, QSizePolicy.Maximum)

        # Label of Combo Box
        label_combo = CustomLabel("File Format: ")
        label_combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
        # ComboBox
        Combo = CustomCombo()
        Combo.addItem("smi format")
        Combo.addItem("txt format")
        Combo.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        FFS2 = QSpacerItem(1000, 10, QSizePolicy.Expanding, QSizePolicy.Expanding)

        # File path
        label_file = CustomLabel("File Path: ")
        label_file.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        label_file.setFixedHeight(30)

        # File Label for path
        label_path = CustomPathLabel("")

        # Button for file path
        Button_path = QPushButton("...")
        Button_path.setFixedWidth(40)
        Button_path.setFixedHeight(40)
        Button_path.clicked.connect(lambda: self.open_file_dialog(self.get_file_filter_44(Combo), label_path))

        # Horizontal layout
        Hlayout = QHBoxLayout()

        # Line Edit
        TextLE = CustomLabel("Name folder: ")
        MyLE = CustomLineEdit()
        validator = QRegularExpressionValidator(QRegularExpression("^[A-Za-z ]*$"))
        MyLE.setValidator(validator)

        # PushButton
        PButton = CustomButton("Predict")
        PButton.clicked.connect(lambda: self.conversion_to_pdbqt_44(label_path, MyLE, Combo, ProgresBar))
        PButton.setFixedWidth(300)

        # Progress Bar
        ProgresBar = QProgressBar()
        ProgresBar.setFixedHeight(50)

        # Convert layout
        Clayout = QHBoxLayout()
        CSL = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        CSR = QSpacerItem(10000, 100, QSizePolicy.Maximum, QSizePolicy.Maximum)
        
        # Bottom Spacer
        spacer_bottom = QSpacerItem(20, 10, QSizePolicy.Minimum, QSizePolicy.Expanding)


        #  ------------- Add Items -------------
        
        # Add title
        title_layout.addItem(TSL)
        title_layout.addWidget(des_label)
        title_layout.addItem(TSR)
        layout.addItem(title_layout)

        # Add Combo
        HCombolayout = QHBoxLayout()
        layout.addSpacing(100)
        HCombolayout.addWidget(label_combo)
        HCombolayout.addWidget(Combo)
        HCombolayout.addItem(FFS2)
        layout.addItem(HCombolayout)
        layout.addSpacing(50)

        # Horizontal label path and button path layout
        layout.addWidget(label_file)
        Hlayout.addWidget(label_path)
        Hlayout.addWidget(Button_path)
        layout.addItem(Hlayout)
        layout.addSpacing(50)

        # Add the line edit
        layout.addWidget(TextLE)
        layout.addWidget(MyLE)
        layout.addSpacing(50)
        
        # Horizontal layout for the button
        Clayout.addItem(CSL)
        Clayout.addWidget(PButton)
        Clayout.addItem(CSR)
        layout.addItem(Clayout)

        # Progress bar and bottom
        layout.addSpacing(20)
        layout.addWidget(ProgresBar)
        layout.addItem(spacer_bottom)

        
        # Set layout and final spacer
        layout.addItem(spacer_bottom)
        self.page_44.setLayout(layout)
    
    def conversion_to_pdbqt_44(self, file_path, des_folder, Combo, bar):

        file_path = file_path.text()
        folder_text = des_folder.text() if des_folder.text() != "" else "BBB predictions"
        Combo = Combo.currentIndex()
        

        if os.path.isfile(file_path):
            
            if Combo == 0:
                try:
                    Con = Consensus_model_for_smi.Conversions()
                    Con.Maximum(file_path)
                    bar.setMaximum(Con.maxim)

                    try:
                        for smile in Con.conversions(file_path, folder_text):
                            bar.setValue(Con.contator)
                    except Exception as e:
                        QMessageBox.critical(self, "Error 4_4-4", str(e))                 
                    
                except Exception as e:
                    QMessageBox.critical(self, "Error 4_4-5", e)

            
            elif Combo == 1:
                try:
                    Con = Consensus_model_for_txt.Conversions()
                    Con.Maximum(file_path)
                    bar.setMaximum(Con.maxim)

                    try:
                        for smile in Con.conversions(file_path, folder_text):
                            bar.setValue(Con.contator)
                    except Exception as e:
                        QMessageBox.critical(self, "Error 4_4-3", str(e))                 
                    
                except Exception as e:
                    QMessageBox.critical(self, "Error 4_4-2", e)
            
        else:
            QMessageBox.critical(self, "Error 4_4-2", "Please fill in all fields")

    def get_file_filter_44(self, combo):
        if combo.currentIndex() == 0:
            return "smi files (*.smi);;All files (*)"
        elif combo.currentIndex() == 1:
            return "txt files (*.txt);;All files (*)"
    
    def setup_page45(self):

        layout = QVBoxLayout()
        label = CustomLabel("Este es un QLabel en la Página 45")
        CB = CustomCheckBox()
        button = QPushButton("Este es un QPushButton en la Página 45") 
        layout.addWidget(label)
        layout.addWidget(CB)
        layout.addWidget(button)
        self.page_45.setLayout(layout)
    
    def setup_page46(self):

        layout = QVBoxLayout()
        label = CustomLabel("Este es un QLabel en la Página 46")
        CB = CustomCheckBox()
        button = QPushButton("Este es un QPushButton en la Página 46") 
        layout.addWidget(label)
        layout.addWidget(CB)
        layout.addWidget(button)
        self.page_46.setLayout(layout)
    
    def setup_page47(self):

        layout = QVBoxLayout()
        label = CustomLabel("Este es un QLabel en la Página 47")
        CB = CustomCheckBox()
        button = QPushButton("Este es un QPushButton en la Página 47") 
        layout.addWidget(label)
        layout.addWidget(CB)
        layout.addWidget(button)
        self.page_47.setLayout(layout)
    
    def setup_page48(self):

        layout = QVBoxLayout()
        label = CustomLabel("Este es un QLabel en la Página 48")
        CB = CustomCheckBox()
        button = QPushButton("Este es un QPushButton en la Página 48") 
        layout.addWidget(label)
        layout.addWidget(CB)
        layout.addWidget(button)
        self.page_48.setLayout(layout)

    # ------------------------- ESTABLISH THE PAGES DIRECTIONS -------------------------"

    def open_file_dialog(self, files, label):

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select file",
            "",
            files # File filter
        )

        # Set the file path
        label.setText(file_path)
    
    def open_folder_dialog(self, label):
        folder_path = QFileDialog.getExistingDirectory(
            self,
            "Select folder"
        )
        
        # Set the folder path
        label.setText(folder_path)

    # ------------------------- ESTABLISH THE PAGES DIRECTIONS -------------------------"

    def page_11_signal(self):
        self.main_stack.setCurrentIndex(0)
        self.page_1_stack.setCurrentIndex(0)
    
    def page_12_signal(self):
        self.main_stack.setCurrentIndex(0)
        self.page_1_stack.setCurrentIndex(1)
    
    def page_21_signal(self):
        self.main_stack.setCurrentIndex(1)
        self.page_2_stack.setCurrentIndex(0)

    def page_31_signal(self):
        self.main_stack.setCurrentIndex(2)
        self.page_3_stack.setCurrentIndex(0)
    
    def page_32_signal(self):
        self.main_stack.setCurrentIndex(2)
        self.page_3_stack.setCurrentIndex(1)
    
    def page_33_signal(self):
        self.main_stack.setCurrentIndex(2)
        self.page_3_stack.setCurrentIndex(2)
    
    def page_41_signal(self):
        self.main_stack.setCurrentIndex(3)
        self.page_4_stack.setCurrentIndex(0)
    
    def page_42_signal(self):
        self.main_stack.setCurrentIndex(3)
        self.page_4_stack.setCurrentIndex(1)
    
    def page_43_signal(self):
        self.main_stack.setCurrentIndex(3)
        self.page_4_stack.setCurrentIndex(2)

    def page_44_signal(self):
        self.main_stack.setCurrentIndex(3)
        self.page_4_stack.setCurrentIndex(3)
    
    def page_45_signal(self):
        self.main_stack.setCurrentIndex(3)
        self.page_4_stack.setCurrentIndex(4)


       
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec()

