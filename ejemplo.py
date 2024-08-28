from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QVBoxLayout, QDialog
from PySide6.QtCore import Qt
from PySide6.QtGui import QCursor

class DropdownWindow(QWidget):
    def __init__(self):
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
        self.resize(185, 150)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Main Window")
        
        # Crear el botón principal
        self.button = QPushButton("Show Dropdown")
        self.button.clicked.connect(self.show_dropdown)

        # Crear la ventana principal
        central_widget = QWidget()
        layout = QVBoxLayout(central_widget)
        layout.addWidget(self.button)
        
        self.setCentralWidget(central_widget)
        
        # Crear la ventana desplegable
        self.dropdown = DropdownWindow()

    def show_dropdown(self):
        # Obtener la posición del botón para mostrar la ventana desplegable debajo de él
        button_rect = self.button.rect()
        button_pos = self.button.mapToGlobal(button_rect.bottomLeft())
        
        self.dropdown.move(button_pos)
        self.dropdown.show()

if __name__ == "__main__":
    app = QApplication([])
    
    main_window = MainWindow()
    main_window.show()
    
    app.exec()
