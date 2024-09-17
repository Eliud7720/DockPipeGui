from PySide6.QtWidgets import QLabel

class CustomTitleLabel(QLabel):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)

        self.setStyleSheet("""
            QLabel {
                color: Black;
                font-family: 'Montserrat';
                font-size: 35px;
                font-weight: 700;
                padding-bottom: 1px;
            }
        """)