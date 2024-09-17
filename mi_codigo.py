from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QFontDatabase

# Crea una instancia de QApplication
app = QApplication([])

# Usa QFontDatabase para obtener las fuentes disponibles
font_database = QFontDatabase()
available_fonts = font_database.families()
print(available_fonts)

# Termina la aplicación
app.quit()
