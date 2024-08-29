from PySide6.QtWidgets import QPushButton, QWidget, QVBoxLayout
from PySide6.QtCore import Qt

class DropdownWindow4(QWidget):
    def __init__(self, main_window=None):
        super().__init__()
        self.setWindowFlags(Qt.FramelessWindowHint | Qt.Popup)
        self.setStyleSheet("background-color: lightgray;")
        
        layout = QVBoxLayout()
        
        # Añadir algunos botones a la ventana desplegable
        button1 = QPushButton("Button 1")
        button2 = QPushButton("Button 2")
        button3 = QPushButton("Button 3")
        
        layout.addWidget(button1)
        layout.addWidget(button2)
        layout.addWidget(button3)
        
        self.setLayout(layout)
        self.resize(200, 150)
        
        self.main_window = main_window

    def focusOutEvent(self, event):
        self.close()

    def mousePressEvent(self, event):
        if not self.rect().contains(event.pos()):
            self.close()
            if self.main_window:
                self.main_window.simulate_button_click_4()
        else:
            super().mousePressEvent(event)
