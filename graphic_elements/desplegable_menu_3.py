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

class DropdownWindow3(QWidget):

    button1Clicked = Signal()
    button2Clicked = Signal()
    button3Clicked = Signal()


    def __init__(self, main_window=None):
        super().__init__()
        self.setWindowFlags(Qt.FramelessWindowHint | Qt.Popup)
        self.setStyleSheet("""
                           DropdownWindow3 {
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
        
        # Añadir algunos botones a la ventana desplegable
        button1 = CustomButton("Redocking")
        button2 = CustomButton("Virtual Screening")
        button3 = CustomButton("Consensus docking")
        
        
        layout.addWidget(button1)
        layout.addWidget(button2)
        layout.addWidget(button3)
        
        self.setLayout(layout)
        self.resize(200, 150)  # Ajusta el tamaño de la ventana desplegable
        
        self.main_window = main_window

        button1.clicked.connect(self.button1_clicked)
        button2.clicked.connect(self.button2_clicked)
        button3.clicked.connect(self.button3_clicked)
    
    def button1_clicked(self):
        self.button1Clicked.emit()
        self.close()
    
    def button2_clicked(self):
        self.button2Clicked.emit()
        self.close()
    
    def button3_clicked(self):
        self.button3Clicked.emit()
        self.close()

    def show(self, button):
        # Obtener la posición del botón y ajustar la ventana desplegable
        button_rect = button.rect()
        button_pos = button.mapToGlobal(button_rect.bottomLeft())
        
        # Ajustar la posición para que la esquina superior izquierda de la ventana
        # esté 10 píxeles por encima de la esquina inferior izquierda del botón
        adjusted_pos = QPoint(button_pos.x(), button_pos.y() - 10)
        
        self.move(adjusted_pos)
        self.resize(button.width(), self.height())  # Establecer el ancho de la ventana desplegable al ancho del botón
        super().show()

    def closeEvent(self, event):
        if self.main_window:
            self.main_window.simulate_button_click_3()
        super().closeEvent(event)

    def focusOutEvent(self, event):
        self.close()

    def mousePressEvent(self, event):
        if not self.rect().contains(event.pos()):
            self.close()
        else:
            super().mousePressEvent(event)
