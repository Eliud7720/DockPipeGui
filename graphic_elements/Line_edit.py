from PySide6.QtWidgets import QLineEdit

class CustomLineEdit(QLineEdit):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setStyleSheet("""
            QLineEdit {
                color: Black;
                font-family: 'Aptos';
                font-size: 15px;
                padding: 10px;
                border: 1px solid black;
            }
        """)