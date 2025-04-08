import sys
import os
import shutil
import subprocess
import time
from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QVBoxLayout,
                             QHBoxLayout, QLabel, QFileDialog, QWidget, QProgressBar,
                             QMessageBox, QTableWidget, QTableWidgetItem, QListWidget,
                             QAbstractItemView, QHeaderView)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
import pandas as pd
import numpy as np
import joblib


class InferenceThread(QThread):
    """Thread for running inference"""
    progress_update = pyqtSignal(int)
    status_update = pyqtSignal(str)
    process_complete = pyqtSignal(pd.DataFrame)  # Will emit results DataFrame
    process_error = pyqtSignal(str)

    def __init__(self, model_path, csv_path):
        super().__init__()
        self.model_path = model_path
        self.csv_path = csv_path
        self.threshold = 0.20
        # Define risk thresholds
        self.low_risk_threshold = 0.032  # 3.2%
        self.high_risk_threshold = 0.25  # 25%

    def run(self):
        try:
            self.status_update.emit("Starting inference processing...")
            self.progress_update.emit(25)

            # Load the CSV data
            df = pd.read_csv(self.csv_path)
            self.status_update.emit("Loaded CSV data...")

            # Keep only the first column and drop the second column
            first_col = df.columns[0]
            id_col = df[[first_col]].copy()

            # Drop the first column for inference (keep all data columns)
            X = df.drop(columns=df.columns[:2])

            self.progress_update.emit(50)
            self.status_update.emit("Loading model...")

            # Load the model
            model = joblib.load(self.model_path)

            self.progress_update.emit(75)
            self.status_update.emit("Running inference...")

            # Get probabilities
            y_prob = model.predict_proba(X)[:, 1]  # Probability of positive class

            # Apply threshold for binary prediction
            y_pred = (y_prob >= self.threshold).astype(int)

            # Create results DataFrame with only the first column
            results = id_col.copy()
            results['prediction'] = y_pred
            results['probability'] = y_prob

            # Add risk category based on probability
            def get_risk_category(prob):
                if prob < self.low_risk_threshold:
                    return "Low Risk"
                elif prob >= self.high_risk_threshold:
                    return "High Risk"
                else:
                    return "Intermediate Risk"

            results['risk_category'] = results['probability'].apply(get_risk_category)

            self.progress_update.emit(100)
            self.status_update.emit("Inference complete!")

            # Save results to CSV in the same directory as the application
            output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inference_results.csv')
            results.to_csv(output_path, index=False)

            # Emit the results dataframe
            self.process_complete.emit(results)

        except Exception as e:
            self.status_update.emit(f"Error during inference: {str(e)}")
            self.process_error.emit(str(e))


class MatlabProcessThread(QThread):
    """Thread for running MATLAB processing without freezing the UI"""
    progress_update = pyqtSignal(int)
    status_update = pyqtSignal(str)
    process_complete = pyqtSignal(str)  # Will emit CSV file path
    process_error = pyqtSignal(str)

    def __init__(self, matlab_exe_path, input_folder_path):
        super().__init__()
        self.matlab_exe_path = matlab_exe_path
        self.input_folder_path = input_folder_path

    def run(self):
        try:
            self.status_update.emit("Starting MATLAB processing for all DICOM files...")
            self.progress_update.emit(10)

            # Change to the directory of the MATLAB executable
            original_dir = os.getcwd()
            matlab_dir = os.path.dirname(self.matlab_exe_path)
            self.status_update.emit(f"Changing to directory: {matlab_dir}")
            os.chdir(matlab_dir)

            # Current directory after change
            current_dir = os.getcwd()
            self.status_update.emit(f"Current directory: {current_dir}")

            self.progress_update.emit(30)

            # Run the MATLAB executable with lower priority
            self.status_update.emit(f"Running MATLAB: {self.matlab_exe_path}")

            try:
                # Use subprocess with explicit arguments instead of shell=True
                # and set lower process priority
                if sys.platform == 'win32':
                    # For Windows, use CREATE_NO_WINDOW to avoid console popup
                    startupinfo = subprocess.STARTUPINFO()
                    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                    startupinfo.wShowWindow = subprocess.SW_HIDE

                    process = subprocess.Popen(
                        [self.matlab_exe_path],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        startupinfo=startupinfo,
                        creationflags=subprocess.BELOW_NORMAL_PRIORITY_CLASS
                    )
                else:
                    # For non-Windows platforms
                    process = subprocess.Popen(
                        [self.matlab_exe_path],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )

                # Wait for process to complete with timeout
                try:
                    stdout, stderr = process.communicate(timeout=300)  # 5 minute timeout
                    self.status_update.emit(f"MATLAB processing finished with return code: {process.returncode}")
                except subprocess.TimeoutExpired:
                    process.kill()
                    stdout, stderr = process.communicate()
                    self.status_update.emit("MATLAB process timed out and was terminated")

            except Exception as matlab_error:
                self.status_update.emit(f"Error executing MATLAB: {str(matlab_error)}")
                raise

            self.progress_update.emit(60)

            # Change back to original directory
            os.chdir(original_dir)
            self.status_update.emit(f"Changed back to original directory: {original_dir}")

            # Wait for potential file operations to complete
            self.status_update.emit("Waiting for file operations to complete...")
            time.sleep(2)
            self.progress_update.emit(80)

            # Check for CSV files
            self.status_update.emit("Checking for CSV files...")
            csv_files = []

            for file in os.listdir(self.input_folder_path):
                if file.lower().endswith('.csv'):
                    file_path = os.path.join(self.input_folder_path, file)
                    csv_files.append((file_path, os.path.getmtime(file_path)))

            if csv_files:
                # Sort by modification time (newest first)
                csv_files.sort(key=lambda x: x[1], reverse=True)
                csv_path = csv_files[0][0]
                self.status_update.emit(f"Found CSV file: {csv_path}")
                self.progress_update.emit(100)
                self.process_complete.emit(csv_path)
            else:
                raise Exception("No CSV file was generated by MATLAB.")

        except Exception as e:
            self.status_update.emit(f"Error in processing: {str(e)}")
            self.process_error.emit(str(e))


class ECGProcessingApp(QMainWindow):
    def __init__(self):
        super().__init__()

        # Get the base directory (where the executable is located)
        if getattr(sys, 'frozen', False):
            # We are running in a PyInstaller bundle
            self.base_dir = os.path.dirname(sys.executable)
        else:
            # We are running in a normal Python environment
            self.base_dir = os.path.dirname(os.path.abspath(__file__))

        # Configuration paths with relative paths
        self.matlab_exe_path = os.path.join(self.base_dir, "ecg_processor.exe")
        self.model_path = os.path.join(self.base_dir, "XGB_model.pkl")
        self.input_folder_path = self.base_dir  # Use the base directory itself for input files
        self.results_df = None

        # Try to clear any existing DCM files on startup
        try:
            self.clear_input_folder()
        except Exception as e:
            print(f"Warning: Could not clear input folder: {e}")

        # Selected file tracking
        self.selected_dcm_files = []  # List to hold multiple files
        self.result_csv_file = None

        # Initialize UI
        self.init_ui()

    def init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle("ECG Analysis Tool - Batch Processing")
        self.setGeometry(300, 300, 900, 700)

        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Status label
        self.status_label = QLabel("Ready. Please upload DICOM ECG files.")
        self.status_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(self.status_label)

        # File selection area
        file_layout = QHBoxLayout()

        # List widget to display selected files
        self.file_list = QListWidget()
        self.file_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        file_layout.addWidget(self.file_list, 3)  # Give it more space

        # Buttons layout
        buttons_layout = QVBoxLayout()

        self.upload_button = QPushButton("Add DICOM Files")
        self.upload_button.clicked.connect(self.upload_dicom_files)
        buttons_layout.addWidget(self.upload_button)

        self.remove_button = QPushButton("Remove Selected")
        self.remove_button.clicked.connect(self.remove_selected_files)
        buttons_layout.addWidget(self.remove_button)

        self.clear_button = QPushButton("Clear All Files")
        self.clear_button.clicked.connect(self.clear_file_list)
        buttons_layout.addWidget(self.clear_button)

        file_layout.addLayout(buttons_layout, 1)
        main_layout.addLayout(file_layout)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        main_layout.addWidget(self.progress_bar)

        # Process buttons
        process_layout = QHBoxLayout()

        self.process_button = QPushButton("Process All Files")
        self.process_button.clicked.connect(self.process_all_files)
        self.process_button.setEnabled(False)
        process_layout.addWidget(self.process_button)

        main_layout.addLayout(process_layout)

        # Results table
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(4)  # ID, Prediction, Probability, Risk Category
        self.results_table.setHorizontalHeaderLabels(["ID", "Prediction", "Probability", "Risk Category"])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        main_layout.addWidget(self.results_table)

        # Configuration info
        config_layout = QVBoxLayout()
        self.matlab_path_label = QLabel(f"MATLAB Executable: {self.matlab_exe_path}")
        self.input_folder_label = QLabel(f"Input Folder: {os.path.abspath(self.input_folder_path)}")
        self.model_path_label = QLabel(
            f"Model Path: {os.path.abspath(self.model_path) if os.path.exists(self.model_path) else 'Not found'}")
        config_layout.addWidget(self.matlab_path_label)
        config_layout.addWidget(self.input_folder_label)
        config_layout.addWidget(self.model_path_label)
        main_layout.addLayout(config_layout)

    def clear_input_folder(self):
        """Remove all DCM files from the input folder"""
        try:
            count = 0
            for file in os.listdir(self.input_folder_path):
                if file.lower().endswith('.dcm'):
                    file_path = os.path.join(self.input_folder_path, file)
                    try:
                        os.remove(file_path)
                        count += 1
                    except PermissionError:
                        print(f"Permission denied when trying to delete {file_path}")
                    except Exception as e:
                        print(f"Could not delete {file_path}: {e}")

            if count > 0:
                print(f"Cleared {count} DCM files from input folder")
        except Exception as e:
            print(f"Error clearing input folder: {e}")

    def upload_dicom_files(self):
        """Handle multiple DICOM file selection"""
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFiles)
        file_paths, _ = file_dialog.getOpenFileNames(
            self, "Select DICOM ECG Files", "", "DICOM Files (*.dcm);;All Files (*)"
        )

        if file_paths:
            for file_path in file_paths:
                # Only add if not already in the list
                if file_path not in self.selected_dcm_files:
                    self.selected_dcm_files.append(file_path)
                    self.file_list.addItem(os.path.basename(file_path))

            # Update UI
            self.status_label.setText(
                f"{len(self.selected_dcm_files)} files selected. Click 'Process All Files' to begin.")
            self.process_button.setEnabled(len(self.selected_dcm_files) > 0)
            self.progress_bar.setValue(0)

    def remove_selected_files(self):
        """Remove selected files from the list"""
        selected_items = self.file_list.selectedItems()
        for item in selected_items:
            file_name = item.text()
            item_row = self.file_list.row(item)

            # Remove from list widget
            self.file_list.takeItem(item_row)

            # Remove from selected files list
            for i, file_path in enumerate(self.selected_dcm_files):
                if os.path.basename(file_path) == file_name:
                    self.selected_dcm_files.pop(i)
                    break

        # Update UI
        self.process_button.setEnabled(len(self.selected_dcm_files) > 0)
        self.status_label.setText(f"{len(self.selected_dcm_files)} files selected.")

    def clear_file_list(self):
        """Clear all files from the list"""
        self.file_list.clear()
        self.selected_dcm_files = []
        self.process_button.setEnabled(False)
        self.status_label.setText("Ready. Please upload DICOM ECG files.")
        self.progress_bar.setValue(0)

    def process_all_files(self):
        """Process all selected files together"""
        if not self.selected_dcm_files:
            QMessageBox.warning(self, "Warning", "Please select at least one DICOM file first.")
            return

        # Check if MATLAB executable exists
        if not os.path.exists(self.matlab_exe_path):
            QMessageBox.critical(self, "Error",
                                 f"MATLAB executable not found at: {self.matlab_exe_path}")
            return

        # First, copy all DICOM files to the input folder
        self.status_label.setText("Copying DICOM files to processing folder...")
        self.progress_bar.setValue(10)

        # Clear any existing DICOM files first
        self.clear_input_folder()

        # Copy all selected files to the input folder
        try:
            for i, file_path in enumerate(self.selected_dcm_files):
                filename = os.path.basename(file_path)
                destination = os.path.join(self.input_folder_path, filename)

                # Copy the file
                shutil.copy2(file_path, destination)

                # Update progress
                progress = int(10 + (i / len(self.selected_dcm_files)) * 20)  # 10% to 30%
                self.progress_bar.setValue(progress)

                # Update status occasionally
                if i % 5 == 0 or i == len(self.selected_dcm_files) - 1:
                    self.status_label.setText(f"Copied {i + 1}/{len(self.selected_dcm_files)} files...")
                    QApplication.processEvents()  # Allow UI to update

            self.status_label.setText(f"All {len(self.selected_dcm_files)} DICOM files copied successfully.")
            self.progress_bar.setValue(30)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error copying files: {str(e)}")
            self.status_label.setText("Error copying files. Process aborted.")
            return

        # Disable buttons during processing
        self.upload_button.setEnabled(False)
        self.remove_button.setEnabled(False)
        self.clear_button.setEnabled(False)
        self.process_button.setEnabled(False)

        # Start MATLAB processing thread
        self.matlab_thread = MatlabProcessThread(self.matlab_exe_path, self.input_folder_path)

        # Connect signals
        self.matlab_thread.progress_update.connect(self.progress_bar.setValue)
        self.matlab_thread.status_update.connect(self.status_label.setText)
        self.matlab_thread.process_complete.connect(self.run_inference)
        self.matlab_thread.process_error.connect(self.processing_error)

        # Start processing
        self.matlab_thread.start()

    def run_inference(self, csv_path):
        """Run inference on the generated CSV file"""
        self.status_label.setText(f"Running inference on all data...")

        # Start inference thread
        self.inference_thread = InferenceThread(self.model_path, csv_path)

        # Connect signals
        self.inference_thread.progress_update.connect(self.progress_bar.setValue)
        self.inference_thread.status_update.connect(self.status_label.setText)
        self.inference_thread.process_complete.connect(self.display_results)
        self.inference_thread.process_error.connect(self.inference_error)

        # Start inference
        self.inference_thread.start()

    def display_results(self, results_df):
        """Display inference results"""
        self.results_df = results_df

        # Update the results table
        row_count = len(results_df)
        self.results_table.setRowCount(row_count)

        # Fill the table with data
        for row_idx, (_, row) in enumerate(results_df.iterrows()):
            # ID column (assuming first column is ID)
            id_item = QTableWidgetItem(str(row.iloc[0]))
            self.results_table.setItem(row_idx, 0, id_item)

            # Prediction column
            prediction_item = QTableWidgetItem(f"{'Positive' if row['prediction'] == 1 else 'Negative'}")
            self.results_table.setItem(row_idx, 1, prediction_item)

            # Probability column
            prob_value = row['probability']
            prob_item = QTableWidgetItem(f"{prob_value:.3f}")
            self.results_table.setItem(row_idx, 2, prob_item)

            # Risk Category column
            risk_item = QTableWidgetItem(row['risk_category'])
            self.results_table.setItem(row_idx, 3, risk_item)

        # Resize table columns and make it visible
        self.results_table.resizeColumnsToContents()

        # Update UI elements
        self.status_label.setText(f"Processing complete. Analyzed {row_count} patients.")
        self.upload_button.setEnabled(True)
        self.remove_button.setEnabled(True)
        self.clear_button.setEnabled(True)
        self.process_button.setEnabled(len(self.selected_dcm_files) > 0)

        # Optional: Display a success message
        QMessageBox.information(self, "Success",
                                f"Analysis completed successfully for {row_count} patients.\n\nResults are displayed in the table and saved to 'inference_results.csv'")

    def processing_error(self, error_message):
        """Handle MATLAB processing errors"""
        # Update UI
        self.result_label = QLabel("Error during processing")
        self.upload_button.setEnabled(True)
        self.remove_button.setEnabled(True)
        self.clear_button.setEnabled(True)
        self.process_button.setEnabled(len(self.selected_dcm_files) > 0)
        self.progress_bar.setValue(0)

        # Display error message
        QMessageBox.critical(self, "Processing Error",
                             f"An error occurred during MATLAB processing:\n{error_message}")

    def inference_error(self, error_message):
        """Handle inference errors"""
        # Update UI
        self.status_label.setText(f"Error during inference")
        self.upload_button.setEnabled(True)
        self.remove_button.setEnabled(True)
        self.clear_button.setEnabled(True)
        self.process_button.setEnabled(len(self.selected_dcm_files) > 0)
        self.progress_bar.setValue(0)

        # Display error message
        QMessageBox.critical(self, "Inference Error",
                             f"An error occurred during inference:\n{error_message}")


if __name__ == "__main__":
    # Set application name and organization for QSettings
    QApplication.setApplicationName("ECG Analysis Tool")
    QApplication.setOrganizationName("Medical Software")

    # Create and show the application
    app = QApplication(sys.argv)
    window = ECGProcessingApp()
    window.show()
    sys.exit(app.exec_())