from PySide6.QtWidgets import QLabel

class CustomLabel(QLabel):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)

        self.setStyleSheet("""
            QLabel {
                color: #424242;
                font-family: 'Aptos';
                font-size: 30px; 
                font-weight: 600;
            }
        """)