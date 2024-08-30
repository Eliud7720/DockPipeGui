from PySide6.QtWidgets import QComboBox

class CustomCombo(QComboBox):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStyleSheet("""
            QComboBox {
                color: black;
                padding: 12px;
                font-size: 20px;
            }

        """)