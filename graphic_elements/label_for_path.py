from PySide6.QtWidgets import QLabel

class CustomPathLabel(QLabel):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)

        self.setStyleSheet("""
            QLabel {
                color: Black;
                font-family: 'Aptos';
                font-size: 15px;
                padding: 10px;
                border: 1px solid black;
            }
        """)