import sys
import subprocess
import os
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QPushButton, QFileDialog, QSizePolicy, QHBoxLayout, QVBoxLayout, QWidget

home_directory = os.path.expanduser("~")
env_directory = "/miniconda3/envs/BiomiX-env/bin/python3"
conda_python_path = "".join([home_directory, env_directory])
print(conda_python_path)

# Get the path of the desired Python executable
desired_python_executable = conda_python_path
print(desired_python_executable)

# Other parts of your script...

def rest_of_the_script():
    class MainWindow(QMainWindow):
        def __init__(self):
            super().__init__()

            self.setWindowTitle("BiomiX")
            icon = QtGui.QIcon()
            icon.addPixmap(
                QtGui.QPixmap("BiomiX_logo3.png"),
                QtGui.QIcon.Mode.Normal,
                QtGui.QIcon.State.Off,
            )
            self.setWindowIcon(icon)
            self.setGeometry(200, 200, 800, 600)

            central_widget = QWidget(self)
            self.setCentralWidget(central_widget)

            self.figure_label = QLabel(self)
            self.figure_label.setAlignment(Qt.AlignCenter)
            self.figure_label.setScaledContents(True)  # Enable pixmap scaling

            self.label = QLabel("Welcome to BiomiX", self)
            font = QtGui.QFont()
            font.setFamily("Arial")
            font.setPointSize(22)
            font.setBold(False)
            font.setWeight(50)
            self.label.setFont(font)
            self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)

            self.label2 = QLabel("Upload the multiomics metadata file\n to start the analysis", self)
            font.setPointSize(16)
            self.label2.setFont(font)
            self.label2.setAlignment(Qt.AlignmentFlag.AlignCenter)

            self.button = QPushButton("Upload", self)
            font.setPointSize(16)
            self.button.setFont(font)
            self.button.clicked.connect(self.open_dialog)

            self.figure_path = "BiomiX_logo3.png"  # Set the path to your figure image

            self.label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self.label2.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self.button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self.figure_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

            # Layout for labels
            label_layout = QVBoxLayout()
            label_layout.addWidget(self.label)
            label_layout.addWidget(self.label2)

            # Layout for button and figure
            bottom_layout = QHBoxLayout()
            bottom_layout.addWidget(self.button)
            bottom_layout.addWidget(self.figure_label)

            # Main layout
            main_layout = QVBoxLayout()
            main_layout.addLayout(label_layout)
            main_layout.addLayout(bottom_layout)

            central_widget.setLayout(main_layout)

            # Upload logo image with a delay
            QTimer.singleShot(100, self.load_figure)

        def resizeEvent(self, event):
            self.adjust_font_size()
            self.load_figure()
            super().resizeEvent(event)

        def adjust_font_size(self):
            # Calcola la dimensione del font basata sulla dimensione della finestra
            width = self.width()
            height = self.height()
            font_size_label = int(height / 15)  # Font size as a fraction of window height
            font_size_label2 = int(height / 20) # Different font size for the second label

            # Set dimension font QLabel
            #Title
            font = self.label.font()
            font.setPointSize(font_size_label)
            self.label.setFont(font)
            #Subtitle
            font = self.label2.font()
            font.setPointSize(font_size_label2)
            self.label2.setFont(font)
            
            
        def load_figure(self):
            pixmap = QPixmap(self.figure_path)
            pixmap = pixmap.scaled(self.figure_label.size(), Qt.AspectRatioMode.KeepAspectRatio)
            self.figure_label.setPixmap(pixmap)
            
#The upper part is shared with the Windows, Linux and Mac OS versions. 
#If there are changes in the interface the changes must be set there.

        def open_dialog(self):
            dialog = QFileDialog()
            dialog.setFileMode(QFileDialog.ExistingFile)
            dialog.setNameFilter("Metadata Files (*.txt *.csv *.tsv)")
            if dialog.exec_():
                file_path = dialog.selectedFiles()[0]
                print(file_path)
                directory = os.path.dirname(file_path)
                self.write_directory(file_path)

                # Execute another script
                script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Biomix_interface_Beta_W7.py")
                python_cmd = [conda_python_path, script_path]
                print(python_cmd)
                subprocess.Popen(python_cmd)

            self.close()

        def write_directory(self, directory):
            script_directory = os.path.dirname(os.path.abspath(__file__))
            output_path = os.path.join(script_directory, "directory.txt")
            with open(output_path, "w") as file:
                file.write(directory)

    if __name__ == '__main__':
        app = QApplication([])
        window = MainWindow()
        window.show()
        app.exec_()

    print("BiomiX interface is loading....")

if __name__ == "__main__":
    # Check if the current Python executable is the desired one
    print(sys.executable != desired_python_executable)
    if sys.executable != desired_python_executable:
        # Use subprocess to run the rest of the script with the desired Python executable
        subprocess.run([desired_python_executable, sys.argv[0]] + sys.argv[1:])
    else:
        rest_of_the_script()

