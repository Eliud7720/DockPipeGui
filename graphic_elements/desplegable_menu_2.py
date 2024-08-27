from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton

class Desplegable_menu(QMainWindow):
    def __init__(self):
        super().__init__()

        # Create the desplegable_menu buttons
        self.button_menu2_1 = QPushButton("Redocking")
        self.button_menu2_2 = QPushButton("Virtual Screening")
        self.button_menu2_3 = QPushButton("Consensus docking")


    
        