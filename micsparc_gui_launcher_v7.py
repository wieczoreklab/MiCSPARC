import sys
import subprocess
import os
from PyQt6.QtWidgets import (
    QApplication, QWidget, QLabel, QPushButton, QLineEdit, QVBoxLayout,
    QHBoxLayout, QFileDialog, QTabWidget, QTextEdit, QSpinBox, QMainWindow
)
from PyQt6.QtCore import QThread, pyqtSignal, Qt


class ScriptRunner(QThread):
    output_signal = pyqtSignal(str)

    def __init__(self, command):
        super().__init__()
        self.command = command

    def run(self):
        process = subprocess.Popen(
            self.command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        for line in process.stdout:
            self.output_signal.emit(line)


class UnifyPhiTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        self.input_file = QLineEdit()
        self.browse_input = QPushButton("Browse")
        self.output_name = QLineEdit()

        self.pfn_input = QSpinBox()
        self.pfn_input.setMinimum(1)

        self.threads_input = QSpinBox()
        self.threads_input.setMinimum(1)
        self.threads_input.setValue(10)

        self.help_button = QPushButton("Help")
        self.run_button = QPushButton("Run Unify Phi")

        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("Input .cs:"))
        file_layout.addWidget(self.input_file)
        file_layout.addWidget(self.browse_input)

        layout.addLayout(file_layout)
        layout.addWidget(QLabel("Output filename:"))
        layout.addWidget(self.output_name)
        layout.addWidget(QLabel("PFN:"))
        layout.addWidget(self.pfn_input)
        layout.addWidget(QLabel("Threads:"))
        layout.addWidget(self.threads_input)
        layout.addWidget(self.help_button)
        layout.addWidget(self.run_button)

        self.setLayout(layout)


class UnifyPsiTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        self.input_file = QLineEdit()
        self.browse_input = QPushButton("Browse")
        self.output_name = QLineEdit()
        self.help_button = QPushButton("Help")
        self.run_button = QPushButton("Run Unify Psi")

        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("Input .cs:"))
        file_layout.addWidget(self.input_file)
        file_layout.addWidget(self.browse_input)

        layout.addLayout(file_layout)
        layout.addWidget(QLabel("Output filename:"))
        layout.addWidget(self.output_name)
        layout.addWidget(self.help_button)
        layout.addWidget(self.run_button)
        self.setLayout(layout)


class ExtrapolateTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        self.input_file = QLineEdit()
        self.browse_input = QPushButton("Browse")
        self.output_name = QLineEdit()
        self.threads_input = QSpinBox()
        self.threads_input.setMinimum(1)
        self.threads_input.setValue(10)
        self.help_button = QPushButton("Help")
        self.run_button = QPushButton("Run Extrapolate Filaments")

        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("Input .cs:"))
        file_layout.addWidget(self.input_file)
        file_layout.addWidget(self.browse_input)

        layout.addLayout(file_layout)
        layout.addWidget(QLabel("Output filename:"))
        layout.addWidget(self.output_name)
        layout.addWidget(QLabel("Threads:"))
        layout.addWidget(self.threads_input)
        layout.addWidget(self.help_button)
        layout.addWidget(self.run_button)
        self.setLayout(layout)


class CreatePFNReferencesTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        self.input_file = QLineEdit()
        self.browse_input = QPushButton("Browse")
        self.recenter = QLineEdit()
        self.apix = QLineEdit()
        self.pfn = QLineEdit()
        self.threads_input = QSpinBox()
        self.threads_input.setMinimum(1)
        self.threads_input.setValue(10)
        self.help_button = QPushButton("Help")
        self.run_button = QPushButton("Run Create PFN References")

        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("Input .cs:"))
        file_layout.addWidget(self.input_file)
        file_layout.addWidget(self.browse_input)

        layout.addLayout(file_layout)
        layout.addWidget(QLabel("New center:"))
        layout.addWidget(self.recenter)
        layout.addWidget(QLabel("Pixel size (A/px):"))
        layout.addWidget(self.apix)
        layout.addWidget(QLabel("Original PFN:"))
        layout.addWidget(self.pfn)
        layout.addWidget(QLabel("Threads:"))
        layout.addWidget(self.threads_input)
        layout.addWidget(self.help_button)
        layout.addWidget(self.run_button)
        self.setLayout(layout)


class AssignPFNsTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        self.input_file = QLineEdit()
        self.browse_input = QPushButton("Browse")
        self.conf_input = QLineEdit("0.7")
        self.output_name = QLineEdit()
        self.help_button = QPushButton("Help")
        self.run_button = QPushButton("Run Assign PFNs")

        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("Input .cs:"))
        file_layout.addWidget(self.input_file)
        file_layout.addWidget(self.browse_input)

        layout.addLayout(file_layout)
        layout.addWidget(QLabel("Confidence Threshold:"))
        layout.addWidget(self.conf_input)
        layout.addWidget(QLabel("Output filename:"))
        layout.addWidget(self.output_name)
        layout.addWidget(self.help_button)
        layout.addWidget(self.run_button)
        self.setLayout(layout)


class SeamSearchTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        self.input_file = QLineEdit()
        self.browse_input = QPushButton("Browse")
        self.recenter = QLineEdit()
        self.recenter_init_px = QLineEdit()
        self.recenter_final_px = QLineEdit()
        self.output_name = QLineEdit()
        self.conf_input = QLineEdit("0.5")
        self.threads_input = QSpinBox()
        self.threads_input.setMinimum(1)
        self.threads_input.setValue(10)
        self.help_button = QPushButton("Help")
        self.run_button = QPushButton("Run Seam Search")

        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("Input .cs:"))
        file_layout.addWidget(self.input_file)
        file_layout.addWidget(self.browse_input)

        layout.addLayout(file_layout)
        layout.addWidget(QLabel("New center (semicolon-separated for multiple PFNs):"))
        layout.addWidget(self.recenter)
        layout.addWidget(QLabel("Initial pixel size (A/px):"))
        layout.addWidget(self.recenter_init_px)
        layout.addWidget(QLabel("Final pixel size (A/px):"))
        layout.addWidget(self.recenter_final_px)
        layout.addWidget(QLabel("Confidence threshold:"))
        layout.addWidget(self.conf_input)
        layout.addWidget(QLabel("Output filename:"))
        layout.addWidget(self.output_name)
        layout.addWidget(QLabel("Threads:"))
        layout.addWidget(self.threads_input)
        layout.addWidget(self.help_button)
        layout.addWidget(self.run_button)
        self.setLayout(layout)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MiCSPARC")
        self.output_dir = os.getcwd()
        self.script_dir = os.getcwd()

        self.tabs = QTabWidget()
        self.extrapolate_tab = ExtrapolateTab()
        self.createpfn_tab = CreatePFNReferencesTab()
        self.assign_tab = AssignPFNsTab()
        self.psi_tab = UnifyPsiTab()
        self.phi_tab = UnifyPhiTab()
        self.seam_tab = SeamSearchTab()

        self.tabs.addTab(self.extrapolate_tab, "Extrapolate Filaments")
        self.tabs.addTab(self.createpfn_tab, "Create PFN References")
        self.tabs.addTab(self.assign_tab, "Assign PFNs")
        self.tabs.addTab(self.psi_tab, "Unify Psi")
        self.tabs.addTab(self.phi_tab, "Unify Phi")
        self.tabs.addTab(self.seam_tab, "Seam Search")

        self.log = QTextEdit()
        self.log.setReadOnly(True)

        info_label = QLabel(
            "<b>🧊 MiCSPARC GUI</b><br>MiCSPARC is a microtubule processing pipeline developed around CryoSPARC to determine structures of decorated and undecorated microtubules.<br>"
            ""
            "<a href='https://github.com/wieczoreklab/MiCSPARC'>View on GitHub</a>"
        )
        info_label.setOpenExternalLinks(True)
        info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        info_label.setStyleSheet("font-size: 13px; margin-bottom: 10px;")

        self.output_dir_field = QLineEdit(self.output_dir)
        self.output_dir_button = QPushButton("Select Output Directory")

        self.script_dir_field = QLineEdit(self.script_dir)
        self.script_dir_button = QPushButton("Select Script Directory")

        out_layout = QHBoxLayout()
        out_layout.addWidget(QLabel("Output directory:"))
        out_layout.addWidget(self.output_dir_field)
        out_layout.addWidget(self.output_dir_button)

        script_layout = QHBoxLayout()
        script_layout.addWidget(QLabel("Script directory:"))
        script_layout.addWidget(self.script_dir_field)
        script_layout.addWidget(self.script_dir_button)

        main_layout = QVBoxLayout()
        main_layout.addWidget(info_label)
        main_layout.addLayout(out_layout)
        main_layout.addLayout(script_layout)
        main_layout.addWidget(self.tabs)
        main_layout.addWidget(QLabel("Log Output:"))
        main_layout.addWidget(self.log)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        self.output_dir_button.clicked.connect(self.select_output_directory)
        self.script_dir_button.clicked.connect(self.select_script_directory)

        for tab in [self.extrapolate_tab, self.createpfn_tab, self.assign_tab, self.psi_tab, self.phi_tab, self.seam_tab]:
            tab.browse_input.clicked.connect(lambda _, t=tab: self.select_file(t.input_file))

        self.extrapolate_tab.run_button.clicked.connect(self.run_extrapolate_script)
        self.createpfn_tab.run_button.clicked.connect(self.run_createpfn_script)
        self.assign_tab.run_button.clicked.connect(self.run_assign_script)
        self.psi_tab.run_button.clicked.connect(self.run_psi_script)
        self.phi_tab.run_button.clicked.connect(self.run_phi_script)
        self.seam_tab.run_button.clicked.connect(self.run_seam_script)

        self.phi_tab.help_button.clicked.connect(self.show_phi_help)
        self.psi_tab.help_button.clicked.connect(self.show_psi_help)
        self.extrapolate_tab.help_button.clicked.connect(self.show_extrapolate_help)
        self.createpfn_tab.help_button.clicked.connect(self.show_createpfn_help)
        self.assign_tab.help_button.clicked.connect(self.show_assign_help)
        self.seam_tab.help_button.clicked.connect(self.show_seam_help)


    def select_file(self, line_edit):
        file, _ = QFileDialog.getOpenFileName(self, "Select Input .cs File", "", "CryoSPARC CS (*.cs)")
        if file:
            line_edit.setText(file)

    def select_output_directory(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir:
            self.output_dir = dir
            self.output_dir_field.setText(dir)

    def select_script_directory(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Script Directory")
        if dir:
            self.script_dir = dir
            self.script_dir_field.setText(dir)

    def build_output_path(self, filename):
        return os.path.join(self.output_dir, filename)

    def build_script_path(self, script_name):
        return os.path.join(self.script_dir, script_name)

    def run_script(self, cmd):
        self.log.append(f"\nRunning: {' '.join(cmd)}\n")
        self.thread = ScriptRunner(cmd)
        self.thread.output_signal.connect(self.log.append)
        self.thread.start()

    def run_phi_script(self):
        cmd = [
            "python", self.build_script_path("csparc_unify_phi.py"),
            "--i", self.phi_tab.input_file.text(),
            "--o", self.build_output_path(self.phi_tab.output_name.text()),
            "--pfn", str(self.phi_tab.pfn_input.value()),
            "--j", str(self.phi_tab.threads_input.value())
        ]
        self.run_script(cmd)

    def run_psi_script(self):
        cmd = [
            "python", self.build_script_path("csparc_unify_psi.py"),
            "--i", self.psi_tab.input_file.text(),
            "--o", self.build_output_path(self.psi_tab.output_name.text())
        ]
        self.run_script(cmd)

    def run_extrapolate_script(self):
        cmd = [
            "python", self.build_script_path("csparc_extrapolate_filaments_hmh.py"),
            "--i", self.extrapolate_tab.input_file.text(),
            "--o", self.build_output_path(self.extrapolate_tab.output_name.text()),
            "--j", str(self.extrapolate_tab.threads_input.value())
        ]
        self.run_script(cmd)

    def run_createpfn_script(self):
        cmd = [
            "python", self.build_script_path("csparc_create_pfn_references.py"),
            "--i", self.createpfn_tab.input_file.text(),
            "--recenter", self.createpfn_tab.recenter.text(),
            "--apix", self.createpfn_tab.apix.text(),
            "--pfn", self.createpfn_tab.pfn.text(),
            "--j", str(self.createpfn_tab.threads_input.value())
        ]
        self.run_script(cmd)

    def run_assign_script(self):
        cmd = [
            "python", self.build_script_path("csparc_assign_pfns.py"),
            "--i", self.assign_tab.input_file.text(),
            "--conf", self.assign_tab.conf_input.text()
        ]
        if self.assign_tab.output_name.text():
            cmd += ["--o", self.build_output_path(self.assign_tab.output_name.text())]
        self.run_script(cmd)

    def run_seam_script(self):
        cmd = [
            "python", self.build_script_path("csparc_seam_search_HMH_09.py"),
            "--i", self.seam_tab.input_file.text(),
            "--recenter", self.seam_tab.recenter.text(),
            "--recenter_init_pxsize", self.seam_tab.recenter_init_px.text(),
            "--recenter_final_pxsize", self.seam_tab.recenter_final_px.text(),
            "--conf", self.seam_tab.conf_input.text(),
            "--j", str(self.seam_tab.threads_input.value())
        ]
        if self.seam_tab.output_name.text():
            cmd += ["--o", self.build_output_path(self.seam_tab.output_name.text())]
        self.run_script(cmd)
    def show_phi_help(self):
        self.run_script(["python", self.build_script_path("csparc_unify_phi.py"), "--help"])

    def show_psi_help(self):
        self.run_script(["python", self.build_script_path("csparc_unify_psi.py"), "--help"])

    def show_extrapolate_help(self):
        self.run_script(["python", self.build_script_path("csparc_extrapolate_filaments_hmh.py"), "--help"])

    def show_createpfn_help(self):
        self.run_script(["python", self.build_script_path("csparc_create_pfn_references.py"), "--help"])

    def show_assign_help(self):
        self.run_script(["python", self.build_script_path("csparc_assign_pfns.py"), "--help"])

    def show_seam_help(self):
        self.run_script(["python", self.build_script_path("csparc_seam_search_HMH_09.py"), "--help"])


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.resize(800, 600)
    window.show()
    sys.exit(app.exec())
