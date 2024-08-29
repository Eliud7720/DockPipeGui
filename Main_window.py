import sys
from PySide6.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QFrame, QLabel, QHBoxLayout, QStackedWidget, QPushButton
from graphic_elements.desplegable_menu_1 import DropdownWindow1
from graphic_elements.desplegable_menu_2 import DropdownWindow2
from graphic_elements.desplegable_menu_3 import DropdownWindow3
from graphic_elements.desplegable_menu_4 import DropdownWindow4
from graphic_elements.label_text import CustomLabel
from graphic_elements.py_toggle import CustomCheckBox
from graphic_elements.up_button_triangle import CustomTriangleButton
from PySide6.QtGui import QIcon, QPixmap
from PySide6.QtCore import Qt


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # ------------------------- ESTABLISH THE MAIN ATTRIBUTES -------------------------"

        # Simulated clicks between buttons
        self.cs_1 = False
        self.cs_2 = False
        self.cs_3 = False
        self.cs_4 = False

        # Simulated clicks between windows
        self.ws_1 = False
        self.ws_2 = False
        self.ws_3 = False
        self.ws_4 = False


        # ------------------------- MAIN WINDOW SETTINGS -------------------------"

        
        self.setWindowTitle("DockPipeGui v. 0.1.1") # Establish the Main Window title
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
                background-color: #4c4c4c;
            }
        """)

        self.main_frame.setStyleSheet("""
         QFrame {
                background-color: white;
            }
        """)
        

        # ------------------------- ESTABLISH THE MAIN WIDGET -------------------------"
        self.setCentralWidget(self.main_frame)
        self.dropdown_1 = DropdownWindow1(self)
        self.dropdown_2 = DropdownWindow2(self)
        self.dropdown_3 = DropdownWindow3(self)
        self.dropdown_4 = DropdownWindow4(self)

        # ------------------------- ESTABLISH THE CONNECTIONS -------------------------"
        
        self.button1.clicked.connect(self.signal_1)
        self.button1.clicked.connect(self.show_dropdown_1)

        self.button2.clicked.connect(self.signal_2)
        self.button2.clicked.connect(self.show_dropdown_2)

        self.button3.clicked.connect(self.signal_3)
        self.button3.clicked.connect(self.show_dropdown_3)

        self.button4.clicked.connect(self.signal_4)
        self.button4.clicked.connect(self.show_dropdown_4)


        # ------------------------- DROPSHOW FUNCTIONS -------------------------"


    def show_dropdown_1(self):
        if not self.ws_1:
            button_rect = self.button1.rect()
            button_pos = self.button1.mapToGlobal(button_rect.bottomLeft())
            self.dropdown_1.move(button_pos)
            self.dropdown_1.show()
    
    
    def show_dropdown_2(self):
        if not self.ws_2:
            button_rect = self.button2.rect()
            button_pos = self.button2.mapToGlobal(button_rect.bottomLeft())
            self.dropdown_2.move(button_pos)
            self.dropdown_2.show()

    def show_dropdown_3(self):
        if not self.ws_3:
            button_rect = self.button3.rect()
            button_pos = self.button3.mapToGlobal(button_rect.bottomLeft())
            self.dropdown_3.move(button_pos)
            self.dropdown_3.show()
    
    def show_dropdown_4(self):
        if not self.ws_4:
            button_rect = self.button4.rect()
            button_pos = self.button4.mapToGlobal(button_rect.bottomLeft())
            self.dropdown_4.move(button_pos)
            self.dropdown_4.show()
    
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

        # ------------------------- ESTABLISH THE RECEPTORS -------------------------"

    def signal_1(self):
        
        # Check simulated clicks

        if self.button2.isChecked() and (self.cs_1 == False):
            self.cs_2 = True
            self.button2.click()

        if self.button3.isChecked() and (self.cs_1 == False):
            self.cs_3 = True
            self.button3.click()

        if self.button4.isChecked() and (self.cs_1 == False):
            self.cs_4 = True
            self.button4.click()
        
        # Return to false the simulated clicks
        self.cs_2 = False
        self.cs_3 = False
        self.cs_4 = False

        self.main_stack.setCurrentIndex(0)
    
    def signal_2(self):
        
        # Check simulated clicks

        if self.button1.isChecked() and (self.cs_2 == False):
            self.cs_1 = True
            self.button1.click()

        if self.button3.isChecked() and (self.cs_2 == False):
            self.cs_3 = True
            self.button3.click()

        if self.button4.isChecked() and (self.cs_2 == False):
            self.cs_4 = True
            self.button4.click()
        
        # Return to false the simulated clicks
        self.cs_1 = False
        self.cs_3 = False
        self.cs_4 = False

        self.main_stack.setCurrentIndex(1)


    def signal_3(self):

        if self.button1.isChecked() and (self.cs_3 == False):
            self.cs_1 = True
            self.button1.click()

        if self.button2.isChecked() and (self.cs_3 == False):
            self.cs_2 = True
            self.button2.click()

        if self.button4.isChecked() and (self.cs_3 == False):
            self.cs_4 = True
            self.button4.click()
        
        # Return to false the simulated clicks
        self.cs_1 = False
        self.cs_2 = False
        self.cs_4 = False

        self.main_stack.setCurrentIndex(2)



    def signal_4(self):
        
        if self.button1.isChecked() and (self.cs_4 == False):
            self.cs_1 = True
            self.button1.click()

        if self.button2.isChecked() and (self.cs_4 == False):
            self.cs_2 = True
            self.button2.click()

        if self.button3.isChecked() and (self.cs_4 == False):
            self.cs_3 = True
            self.button3.click()
        
        # Return to false the simulated clicks
        self.cs_1 = False
        self.cs_2 = False
        self.cs_3 = False

        self.main_stack.setCurrentIndex(3)
    


        # ------------------------- ESTABLISH THE PAGES FUNCTIONS -------------------------"

    def setup_page1(self):

        layout = QVBoxLayout()

        label = CustomLabel("Este es un QLabel en la Página 1")
        label2 = CustomLabel("Este es un QLabel en la Página 1, label 2")
        CB = CustomCheckBox()
        button = QPushButton("Este es un QPushButton en la Página 1")

        layout.addWidget(label)
        layout.addWidget(label2)
        layout.addWidget(CB)
        layout.addWidget(button)


        self.page1.setLayout(layout)
    
    def setup_page2(self):
        layout = QVBoxLayout()
        layout.addWidget(QPushButton("Botón en Página 2"))
        self.page2.setLayout(layout)

    def setup_page3(self):
        layout = QVBoxLayout()
        layout.addWidget(QPushButton("Botón en Página 3"))
        self.page3.setLayout(layout)

    def setup_page4(self):
        layout = QVBoxLayout()
        layout.addWidget(QPushButton("Botón en Página 4"))
        self.page4.setLayout(layout)


        
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec()

