import sys
from PySide6.QtWidgets import QDialog, QVBoxLayout, QCheckBox, QDialogButtonBox, QLabel


class CheckableOptionsDialog(QDialog):
    def __init__(self, options, protein, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select the chains")
        
        # Dialog layout
        layout = QVBoxLayout(self)
        
        # Create a label
        self.label = QLabel(f"Please, select the desired chains of {protein}:", self)
        layout.addWidget(self.label)

        # Create a dictionary to store options and checkboxes
        self.checkboxes = {}
        
        # Create checkboxes for each option
        for option in options:
            checkbox = QCheckBox(option, self)
            layout.addWidget(checkbox)
            self.checkboxes[option] = checkbox
        
        # OK button
        buttons = QDialogButtonBox(QDialogButtonBox.Ok, self)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)

    def get_checked_options(self):
        return [option for option, checkbox in self.checkboxes.items() if checkbox.isChecked()]