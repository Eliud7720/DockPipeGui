from PySide6.QtWidgets import QPushButton, QWidget, QVBoxLayout
from PySide6.QtCore import Qt, QPoint
from PySide6.QtCore import Signal

class CustomButton(QPushButton):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.setStyleSheet("""
            QPushButton {
                color: #4c4c4c;
                font-family: 'Aptos';
                font-size: 25px; 
                background-color: white;
                border-radius: 13px;
                border: 1px solid black;
            }
                           
            QPushButton:hover {
                color: #4c4c4c;
                background-color: #c7c7c7; 
            }
                           
            QPushButton:pressed {
                color: #4c4c4c;
                background-color: #a1a1a1;
            }
        """)

class DropdownWindow2(QWidget):

    button1Clicked = Signal()

    def __init__(self, main_window=None):
        super().__init__()
        self.setWindowFlags(Qt.FramelessWindowHint | Qt.Popup)
        self.setStyleSheet("""
            DropdownWindow2 {
                background: qlineargradient(
                    x1: 0, y1: 0,
                    x2: 0, y2: 1,
                    stop: 0 #a1a1a1, stop: 1 #ffffff
                );
                border: 10px solid black;
            }

            QPushButton {
                color: black;
                background-color: white;
            }
            """)
        
        
        layout = QVBoxLayout()
        
        button1 = CustomButton("Openbabel")
        
        
        layout.addWidget(button1)
        
        self.setLayout(layout)
        self.resize(200, 150)
        
        self.main_window = main_window
        button1.clicked.connect(self.button1_clicked)
    
    def button1_clicked(self):
        self.button1Clicked.emit()
        self.close()
    

    def button2_clicked(self):
        self.button2Clicked.emit()
        self.close()


    def show(self, button):
        button_rect = button.rect()
        button_pos = button.mapToGlobal(button_rect.bottomLeft())
        adjusted_pos = QPoint(button_pos.x(), button_pos.y() - 10)
        
        self.move(adjusted_pos)
        self.resize(button.width(), self.height())
        super().show()

    def closeEvent(self, event):
        if self.main_window:
            self.main_window.simulate_button_click_2()
        super().closeEvent(event)

    def focusOutEvent(self, event):
        self.close()

    def mousePressEvent(self, event):
        if not self.rect().contains(event.pos()):
            self.close()
        else:
            super().mousePressEvent(event)
