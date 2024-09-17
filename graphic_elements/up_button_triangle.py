
from PySide6.QtWidgets import QPushButton
from PySide6.QtGui import QPainter, QPolygon, QBrush, QColor, QTransform
from PySide6.QtCore import QPoint, Qt, QPropertyAnimation, QObject, Property, QEasingCurve

class RotatingTriangle(QObject):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._angle = 0

    def get_angle(self):
        return self._angle

    def set_angle(self, value):
        if value != self._angle:
            self._angle = value
            self.parent().update()

    angle = Property(float, get_angle, set_angle)

class CustomTriangleButton(QPushButton):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.is_pressed = False

        self.setCheckable(True)

        self.rotating_triangle = RotatingTriangle(self)

        self.animation = QPropertyAnimation(self.rotating_triangle, b"angle")
        self.animation.setDuration(70)
        self.animation.setEasingCurve(QEasingCurve.Linear)

        self.setStyleSheet("""
            QPushButton {
                color: #4c4c4c;
                padding-top: 12px;
                padding-bottom: 12px;
                font-family: 'Arial';
                font-size: 35px; 
                background-color: white;
                border-radius: 13px;
                font-weight: 600;
            }
                           
            QPushButton:hover {
                color: #4c4c4c;
                background-color: #c7c7c7; 
            }
                           
            QPushButton:pressed {
                color: #4c4c4c;
                background-color: #a1a1a1;
            }
            
            QPushButton:checked {
                color: #4c4c4c;
                background-color: #a1a1a1;
            }
        """)

        self.clicked.connect(self.on_clicked)  # Conectar la señal clicked al nuevo método

    def paintEvent(self, event):
        super().paintEvent(event)

        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        triangle_size = 15
        x_pos = self.width() - triangle_size - 10
        y_pos = (self.height() - triangle_size) / 2

        points = [
            QPoint(x_pos, y_pos + triangle_size),  
            QPoint(x_pos + triangle_size, y_pos + triangle_size),  
            QPoint(x_pos + triangle_size / 2, y_pos)  
        ]
        triangle = QPolygon(points)

        brush = QBrush(QColor("#4c4c4c"))
        painter.setBrush(brush)
        painter.setPen(Qt.NoPen)

        transform = QTransform()
        transform.translate(x_pos + triangle_size / 2, y_pos + triangle_size / 2)
        transform.rotate(self.rotating_triangle.angle)
        transform.translate(-(x_pos + triangle_size / 2), -(y_pos + triangle_size / 2))
        painter.setTransform(transform, True)

        painter.drawPolygon(triangle)

    def on_clicked(self):
        if self.isChecked():
            self.animation.setStartValue(self.rotating_triangle.angle)
            self.animation.setEndValue(180)
        else:
            self.animation.setStartValue(self.rotating_triangle.angle)
            self.animation.setEndValue(0)

        self.animation.start()

    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
