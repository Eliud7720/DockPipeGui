from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton, QWidget, QVBoxLayout
from PySide6.QtCore import Qt
from PySide6.QtGui import QCursor

class DropdownWindow(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowFlags(Qt.FramelessWindowHint | Qt.Popup)
        self.setStyleSheet("background-color: lightgray;")
        
        layout = QVBoxLayout()

        button1 = QPushButton("Button 1")
        button2 = QPushButton("Button 2")
        button3 = QPushButton("Button 3")

        layout.addWidget(button1)
        layout.addWidget(button2)
        layout.addWidget(button3)
        
        self.setLayout(layout)
        self.resize(185, 150)
        self.main_window = None  # Aquí guardaremos una referencia al MainWindow para invocar el clic en self.button

    def mousePressEvent(self, event):
        # Si se hace clic fuera del área de la ventana desplegable, se cierra y se simula un clic en el botón de la ventana principal
        if not self.rect().contains(event.pos()):
            self.close()
            if self.main_window:
                self.main_window.simulate_button_click()
        else:
            super().mousePressEvent(event)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Main Window")
        
        self.button = QPushButton("Show Dropdown")
        self.button.clicked.connect(self.show_dropdown)

        central_widget = QWidget()
        layout = QVBoxLayout(central_widget)
        layout.addWidget(self.button)
        
        self.setCentralWidget(central_widget)
        
        self.dropdown = DropdownWindow()
        self.dropdown.main_window = self  # Referencia a la ventana principal

    def show_dropdown(self):
        button_rect = self.button.rect()
        button_pos = self.button.mapToGlobal(button_rect.bottomLeft())
        
        self.dropdown.move(button_pos)
        self.dropdown.show()

    def simulate_button_click(self):
        # Simular un clic en el botón
        self.button.click()

if __name__ == "__main__":
    app = QApplication([])
    
    main_window = MainWindow()
    main_window.show()
    
    app.exec()
