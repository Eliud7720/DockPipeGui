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
                background-color: #00C0EC;
                border-radius: 13px
            }
                           
            QPushButton:hover {
                color: white;
                background-color: #00B0D8; 
            }
                           
            QPushButton:pressed {
                color: white;
                background-color: #0085A5;
            }
        """)


