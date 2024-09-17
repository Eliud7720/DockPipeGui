from PySide6.QtWidgets import QPushButton

class CustomButton(QPushButton):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.setStyleSheet("""
            QPushButton {
                color: white;
                padding: 12px;
                font-family: 'Aptos';
                font-size: 35px; 
                background-color: #207c7a;
                border-radius: 13px;
                font-weight: 600;
            }
                           
            QPushButton:hover {
                color: white;
                background-color: #106a69; 
            }
                           
            QPushButton:pressed {
                color: white;
                background-color: #005958;
            }
        """)


