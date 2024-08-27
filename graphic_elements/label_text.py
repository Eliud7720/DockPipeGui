from PySide6.QtWidgets import QLabel

class CustomLabel(QLabel):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)

        self.setStyleSheet("""
            QLabel {
                color: Black;
                font-family: 'Aptos';
                font-size: 20px; 
            }
        """)