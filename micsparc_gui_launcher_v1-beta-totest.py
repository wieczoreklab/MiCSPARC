import sys
import subprocess
import os
from PyQt6.QtWidgets import (
    QApplication, QWidget, QLabel, QPushButton, QLineEdit, QVBoxLayout,
    QHBoxLayout, QFileDialog, QTabWidget, QTextEdit, QSpinBox, QMainWindow
)
from PyQt6.QtCore import QThread, pyqtSignal, Qt
from PyQt6.QtCore import QFileInfo
from PyQt6.QtWidgets import QFormLayout
from PyQt6.QtWidgets import QMessageBox

class ScriptRunnerWithCWD(QThread):
    output_signal = pyqtSignal(str)

    def __init__(self, command, cwd=None):
        super().__init__()
        self.command = command
        self.cwd = cwd

    def run(self):
        process = subprocess.Popen(
            self.command,
            cwd=self.cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
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

class ImportResultGroupTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        self.project_input = QLineEdit()
        self.workspace_input = QLineEdit()
        self.csg_file_input = QLineEdit()
        self.browse_button = QPushButton("Browse .csg File")
        self.run_button = QPushButton("Run Import")
        self.help_button = QPushButton("Help")

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)

        form_layout = QFormLayout()
        form_layout.addRow("Project UID:", self.project_input)
        form_layout.addRow("Workspace UID:", self.workspace_input)

        file_layout = QHBoxLayout()
        file_layout.addWidget(self.csg_file_input)
        file_layout.addWidget(self.browse_button)

        layout.addLayout(form_layout)
        layout.addWidget(QLabel(".csg file path:"))
        layout.addLayout(file_layout)
        layout.addWidget(self.run_button)
        layout.addWidget(QLabel("Log:"))
        layout.addWidget(self.log_box)
        layout.addWidget(self.help_button)

        self.setLayout(layout)

        self.browse_button.clicked.connect(self.browse_file)
        self.run_button.clicked.connect(self.run_import)
        self.help_button.clicked.connect(self.show_help)

    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select .csg File", "", "CryoSPARC Result Group (*.csg)")
        if file_path:
            self.csg_file_input.setText(file_path)

    def run_import(self):
        from cryosparc.tools import CryoSPARC
        import os
        import yaml

        try:
            self.log_box.append("🔄 Reading config and connecting to CryoSPARC...")
            script_dir = os.path.dirname(os.path.realpath(__file__))
            config_path = os.path.join(script_dir, "cs_config.yml")
            with open(config_path, "r") as f:
                cs_conf = yaml.safe_load(f)

            cryo = CryoSPARC(
                license=cs_conf["license"],
                host=cs_conf["host"],
                base_port=cs_conf["base_port"],
                email=cs_conf["email"],
                password=cs_conf["password"]
            )

            assert cryo.test_connection()
            self.log_box.append("✅ Connected to CryoSPARC.")

            project_uid = self.project_input.text().strip()
            workspace_uid = self.workspace_input.text().strip()

            self.log_box.append(f"📁 Using project: {project_uid}")
            self.log_box.append(f"📁 Using workspace: {workspace_uid}")

            project = cryo.find_project(project_uid)

            job = project.create_job(
                workspace_uid=workspace_uid,
                type="import_result_group",
                title=f"Imported result group from {os.path.basename(self.csg_file_input.text())}",
                params={"blob_path": self.csg_file_input.text()}
            )
            job.queue()

            self.log_box.append(f"✅ Import job {str(job.uid)} submitted.")
        except Exception as e:
            self.log_box.append(f"❌ Error: {str(e)}")

    def show_help(self):
        self.log_box.append(
            "Ensure that the script directory contains a valid cs_config.yml file with the correct CryoSPARC instance credentials. "
            "Also, make sure the installed version of cryosparc-tools matches your CryoSPARC major and minor version. "
            "For example, if your CryoSPARC version is v4.7.1, use the latest v4.7.x release of cryosparc-tools. "
            "Then, browse for your MiCSPARC result file (*.csg) and specify the CryoSPARC project and workspace where you'd like to import the result group."
        )

class AnalyseResults(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        self.pdf_file_field = QLineEdit()
        self.browse_button = QPushButton("Browse .PDF File")
        self.analyse_button = QPushButton("Analyse Results")
        self.help_button = QPushButton("Help")

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)

        file_layout = QHBoxLayout()
        file_layout.addWidget(self.pdf_file_field)
        file_layout.addWidget(self.browse_button)
        file_layout.addWidget(self.analyse_button)

        layout.addWidget(QLabel("Results file path:"))
        layout.addLayout(file_layout)
        layout.addWidget(self.help_button)
        layout.addWidget(QLabel("Log:"))
        layout.addWidget(self.log_box)

        self.setLayout(layout)

        self.browse_button.clicked.connect(self.browse_file)
        self.analyse_button.clicked.connect(self.analyse_results)
        self.help_button.clicked.connect(self.show_help)

    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select .PDF File", "", "MiCSPARC Result (*.pdf)")
        if file_path:
            self.pdf_file_field.setText(file_path)

    def analyse_results(self):
        import os
        from PyQt6.QtGui import QDesktopServices
        from PyQt6.QtCore import QUrl

        file_path = self.pdf_file_field.text()

        if not file_path:
            self.log_box.append(
                "No MiCSPARC Result (*.pdf) found - Please browse to the result file from Assign PFNs, Unify Phi, Unify Psi, or Seam Search."
            )
            return

        if not os.path.isfile(file_path):
            self.log_box.append(f"PDF file not found: {file_path}")
            return

        self.log_box.append(f"Opening PDF: {file_path}")
        QDesktopServices.openUrl(QUrl.fromLocalFile(file_path))

    def show_help(self):
        self.log_box.append(
            "Browse to a MiCSPARC result file (*.pdf) generated by one of the processing steps: Assign PFNs, Unify Phi, Unify Psi, or Seam Search. These files summarize the outcome of each respective analysis and are located in the corresponding output directories."
        )

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
        self.tabs.addTab(self.createpfn_tab, "Create References")
        self.tabs.addTab(self.assign_tab, "Assign PFNs")
        self.tabs.addTab(self.psi_tab, "Unify Psi")
        self.tabs.addTab(self.phi_tab, "Unify Phi")
        self.tabs.addTab(self.seam_tab, "Seam Search")

        # Add extra tabs
        self.import_tab = ImportResultGroupTab()
        self.tabs.addTab(self.import_tab, "Import Result Group into CryoSPARC")
        self.import_tab = AnalyseResults()
        self.tabs.addTab(self.import_tab, "Analyse Results")

        self.log = QTextEdit()
        self.log.setReadOnly(True)

        info_label = QLabel(
            "<b>🧊 MiCSPARC GUI 🧊</b><br><br><br>"
            "<a>MiCSPARC is a Microtubule processing pipeline developed around CryoSPARC to determine structures of decorated and undecorated microtubules.</a><br>"
            "<a href='https://github.com/wieczoreklab/MiCSPARC'>View on GitHub</a><br>"
            "<a href='https://github.com/wieczoreklab/MiCSPARC'>If MiCSPARC contributed to your work, please cite it in your publication.</a><br>"
        )

        info_label.setOpenExternalLinks(True)
        info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        info_label.setStyleSheet("font-size: 16px; margin-bottom: 9px;")

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
            info = QFileInfo(file)
            real_path = info.symLinkTarget() if info.isSymLink() else file
            line_edit.setText(real_path)

    def select_output_directory(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir:
            info = QFileInfo(dir)
            real_dir = info.symLinkTarget() if info.isSymLink() else dir
            self.output_dir = real_dir
            self.output_dir_field.setText(real_dir)

    def select_script_directory(self):
        dir = QFileDialog.getExistingDirectory(self, "Select Script Directory")
        if dir:
            info = QFileInfo(dir)
            real_dir = info.symLinkTarget() if info.isSymLink() else dir
            self.script_dir = real_dir
            self.script_dir_field.setText(real_dir)

    def build_output_path(self, filename):
        return os.path.join(self.output_dir, filename)

    def build_script_path(self, script_name):
        return os.path.join(self.script_dir, script_name)

    def run_script(self, cmd, cwd=None):
        self.log.append(f"\nRunning: {' '.join(cmd)} in {cwd if cwd else os.getcwd()}\n")
        self.thread = ScriptRunnerWithCWD(cmd, cwd=cwd)
        self.thread.output_signal.connect(self.log.append)
        self.thread.start()


    def run_phi_script(self):
        cmd = [
            "python", self.build_script_path("csparc_unify_phi_GUI.py"),
            "--i", os.path.basename(self.phi_tab.input_file.text()),
            "--o", os.path.basename(self.phi_tab.output_name.text()),
            "--pfn", str(self.phi_tab.pfn_input.value()),
            "--j", str(self.phi_tab.threads_input.value())
        ]
        self.run_script(cmd, cwd=self.output_dir)

    def run_psi_script(self):
        cmd = [
            "python", self.build_script_path("csparc_unify_psi_GUI.py"),
            "--i", os.path.basename(self.psi_tab.input_file.text()),
            "--o", os.path.basename(self.psi_tab.output_name.text())
        ]
        self.run_script(cmd, cwd=self.output_dir)

    def run_extrapolate_script(self):
        cmd = [
            "python", self.build_script_path("csparc_extrapolate_filaments_GUI.py"),
            "--i", os.path.basename(self.extrapolate_tab.input_file.text()),
            "--o", os.path.basename(self.extrapolate_tab.output_name.text()),
            "--j", str(self.extrapolate_tab.threads_input.value())
        ]
        self.run_script(cmd, cwd=self.output_dir)

    def run_createpfn_script(self):
        cmd = [
            "python", self.build_script_path("csparc_create_pfn_references_GUI.py"),
            "--i", os.path.basename(self.createpfn_tab.input_file.text()),
            "--recenter", self.createpfn_tab.recenter.text(),
            "--apix", self.createpfn_tab.apix.text(),
            "--pfn", self.createpfn_tab.pfn.text(),
            "--j", str(self.createpfn_tab.threads_input.value())
        ]
        self.run_script(cmd, cwd=self.output_dir)

    def run_assign_script(self):
        cmd = [
            "python", self.build_script_path("csparc_assign_pfns_GUI.py"),
            "--i", os.path.basename(self.assign_tab.input_file.text()),
            "--conf", self.assign_tab.conf_input.text()
        ]
        if self.assign_tab.output_name.text():
            cmd += ["--o", os.path.basename(self.assign_tab.output_name.text())]
        self.run_script(cmd, cwd=self.output_dir)

    def run_seam_script(self):
        cmd = [
            "python", self.build_script_path("csparc_seam_search_GUI.py"),
            "--i", os.path.basename(self.seam_tab.input_file.text()),
            "--recenter", self.seam_tab.recenter.text(),
            "--recenter_init_pxsize", self.seam_tab.recenter_init_px.text(),
            "--recenter_final_pxsize", self.seam_tab.recenter_final_px.text(),
            "--conf", self.seam_tab.conf_input.text(),
            "--j", str(self.seam_tab.threads_input.value())
        ]
        if self.seam_tab.output_name.text():
            cmd += ["--o", os.path.basename(self.seam_tab.output_name.text())]
        self.run_script(cmd, cwd=self.output_dir)

    def show_phi_help(self):
        self.run_script(["python", self.build_script_path("csparc_unify_phi_GUI.py"), "--help"], cwd=self.output_dir)

    def show_psi_help(self):
        self.run_script(["python", self.build_script_path("csparc_unify_psi_GUI.py"), "--help"], cwd=self.output_dir)

    def show_extrapolate_help(self):
        self.run_script(["python", self.build_script_path("csparc_extrapolate_filaments_GUI.py"), "--help"], cwd=self.output_dir)

    def show_createpfn_help(self):
        self.run_script(["python", self.build_script_path("csparc_create_pfn_references_GUI.py"), "--help"], cwd=self.output_dir)

    def show_assign_help(self):
        self.run_script(["python", self.build_script_path("csparc_assign_pfns_GUI.py"), "--help"], cwd=self.output_dir)

    def show_seam_help(self):
        self.run_script(["python", self.build_script_path("csparc_seam_search_GUI.py"), "--help"], cwd=self.output_dir)

def cs_import_result_group(cryo, project_uid, workspace_uid, csg_path, title):
    project = cryo.get_project(project_uid)
    job = project.create_job(
        workspace_uid=workspace_uid,
        type="import_result_group",
        title=title,
        params={"blob_path": csg_path}
    )
    job.queue()
    return job

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.resize(1080, 1080)
    window.show()
    sys.exit(app.exec())
