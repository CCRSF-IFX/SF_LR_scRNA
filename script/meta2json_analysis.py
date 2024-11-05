#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:57:49 2024
@author: ccrsfifx
"""

import os
import json
import tarfile
import argparse
import subprocess
import pandas as pd
import logging
from datetime import datetime, date

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Global variables
fmetadata = None
fanalysis = None
fproject = None
in_path = os.getcwd()

parser = argparse.ArgumentParser(description="""Process Analysis_Combined directory and create sample folders with metadata.""")
parser.add_argument("-m", "--metadata", metavar="metadata.txt", dest="fmetadata", type=str, required=True, help="Input metadata.txt")
parser.add_argument("-af", "--analysis_folder", metavar="analysis_folder", dest="analysisFolder", action="store", type=str, required=True, help="Path to Snakemake root directory containing .library.csv files.")
parser.add_argument("-p", "--projectName", metavar="project_name", dest="fprojectName", type=str, required=True, help="Project name")
args = parser.parse_args()

fmetadata = args.fmetadata
fanalysis = args.analysisFolder.rstrip("/")
fproject = args.fprojectName

class Sample:
    def __init__(self, attributes, values):
        self.attribute2value = {attr: val if val else "Unspecified" for attr, val in zip(attributes, values)}
        self.flag = 1

    def add_attribute(self, attribute, value):
        self.attribute2value[attribute] = value

    def get_attribute(self, attribute):
        return self.attribute2value.get(attribute, "Unspecified")

def load_metadata():
    dtype = {'Sample name': str}
    data = pd.read_csv(fmetadata, header=0, delimiter=',', dtype=dtype)
    attributes = ['PIName', 'PIfirstName', 'PIlastName', 'LabContact', 'BAContact', 'Project_date', 'Platform', 'Application']
    values = []
    
    PIName = data['Principale Investigator'][0]
    PIfirstName = PIName.split(" ")[0]
    PIlastName = PIName.split(" ")[-1]
    values.extend([PIName, PIfirstName, PIlastName])
    values.append(data['Lab contact'][0])
    values.append(data['Bioinformatis contact'][0])

    date_string = data['Run start date '][0]
    date_formats = ['%m/%d/%y', '%m/%d/%Y', '%Y-%m-%d']
    Project_date = None
    for date_format in date_formats:
        try:
            Project_date = datetime.strptime(date_string.split()[0], date_format).strftime('%Y-%m-%d')
            break
        except ValueError:
            continue

    if Project_date is None:
        logging.error(f"Date format for '{date_string}' is not recognized.")
        raise ValueError("Unrecognized date format.")

    values.append(Project_date)
    values.append(data['FC product code'][0])
    values.append(data['Application'][0])

    sample = Sample(attributes, values)
    logging.info(f"Loaded sample metadata with Project date: {sample.get_attribute('Project_date')}")
    return sample

def read_sample_files():
    sample_files = [csv_file for csv_file in os.listdir(fanalysis) if csv_file.endswith(".library.csv")]
    samplerun = []
    for library_file in sample_files:
        with open(os.path.join(fanalysis, library_file), 'r') as file:
            for line in file.readlines():
                if not line.startswith("Sample_ID,Sample,run,path"):
                    sample_now, _, run_now, path = line.rstrip().split(",")[0:4]
                    samplerun.append(f"{sample_now}|{run_now}")
    return list(set(samplerun))

def create_pilab_json(pi_path, pi_name):
    metadata = {
        "metadataEntries": [
            {"attribute": "collection_type", "value": "PI_Lab"},
            {"attribute": "data_owner_designee", "value": pi_name},
            {"attribute": "data_generator", "value": "Zhao, Yongmei"},
            {"attribute": "data_owner", "value": pi_name}
        ]
    }
    with open(f"{pi_path}.metadata.json", "w") as out:
        json.dump(metadata, out, indent=4)
    logging.info(f"Created PI lab metadata.json for {pi_name} at {pi_path}")

def create_project_json(project_path, project_name, LabContact, BAContact, Project_date):
    nas_no = project_name.split("_")[1]
    today = date.today().isoformat()
    metadata = {
        "metadataEntries": [
            {"attribute": "collection_type", "value": "Project"},
            {"attribute": "project_id_CSAS_NAS", "value": nas_no},
            {"attribute": "project_name", "value": project_name},
            {"attribute": "contact_name", "value": LabContact},
            {"attribute": "bioinformatics_contact", "value": BAContact},
            {"attribute": "project_start_date", "value": Project_date},
            {"attribute": "project_completed_date", "value": today},
            {"attribute": "grant_funding_agent", "value": "NIH"}
        ]
    }
    with open(f"{project_path}.metadata.json", "w") as out:
        json.dump(metadata, out, indent=4)
    logging.info(f"Created project metadata.json for {project_name} at {project_path}")

def create_sample_json(sample_path, sample_name):
    metadata = {
        "metadataEntries": [
            {"attribute": "collection_type", "value": "Sample"},
            {"attribute": "sample_id", "value": sample_name},
            {"attribute": "sample_name", "value": sample_name},
            {"attribute": "initial_sample_concentration_ngul", "value": "Unspecified"},
            {"attribute": "initial_sample_volume_ul", "value": "Unspecified"},
            {"attribute": "library_id", "value": "Unspecified"},
            {"attribute": "library_lot", "value": "Unspecified"},
            {"attribute": "sfqc_library_concentration_nM", "value": "Unspecified"},
            {"attribute": "sfqc_library_size", "value": "Unspecified"},
            {"attribute": "source_id", "value": "Unspecified"},
            {"attribute": "source_name", "value": "Unspecified"},
            {"attribute": "source_organism", "value": "Unspecified"},
            {"attribute": "source_provider", "value": "Unspecified"}
        ]
    }
    with open(f"{sample_path}.metadata.json", "w") as out:
        json.dump(metadata, out, indent=4)
    logging.info(f"Created sample metadata.json for {sample_name} at {sample_path}")

def create_flowcell_json(flowcell_path, flowcell_id, instrument, run_date, runid, Application):
    metadata = {
        "metadataEntries": [
            {"attribute": "collection_type", "value": "Flowcell"},
            {"attribute": "flowcell_id", "value": flowcell_id},
            {"attribute": "run_name", "value": runid},
            {"attribute": "run_date", "value": run_date},
            {"attribute": "sequencing_platform", "value": instrument},
            {"attribute": "sequencing_application_type", "value": Application},
            {"attribute": "read_length", "value": "Unspecified"},
            {"attribute": "pooling", "value": "Unspecified"}
        ]
    }
    with open(f"{flowcell_path}.metadata.json", "w") as out:
        json.dump(metadata, out, indent=4)
    logging.info(f"Created flowcell metadata.json at {flowcell_path}")

def create_object_json(sample_path, file_name, samplename, file_type="TAR", compression="Compressed", software_tool="N/A"):
    metadata = {
        "metadataEntries": [
            {"attribute": "object_name", "value": file_name},
            {"attribute": "file_type", "value": file_type},
            {"attribute": "reference_genome", "value": "N/A"},
            {"attribute": "reference_annotation", "value": "N/A"},
            {"attribute": "software_tool", "value": software_tool},
            {"attribute": "data_compression_status", "value": compression},
            {"attribute": "created_by", "value": samplename},
            {"attribute": "date_created", "value": datetime.now().isoformat()}
        ]
    }
    metadata_file_path = os.path.join(sample_path, f"{file_name}.metadata.json")
    with open(metadata_file_path, 'w') as json_file:
        json.dump(metadata, json_file, indent=4)
    logging.info(f"Created metadata.json for {file_name} at {metadata_file_path}")

def config_working_directory(project_name, sample):
    path_hpc_dme_clu = '/mnt/ccrsf-ifx/Software/tools/HPC_DME_APIs'

    with open("dm_register_analysis.sh", "w") as out:
        out.write("#!/bin/bash\n")
        out.write("#SBATCH --partition=norm\n")
        out.write("#SBATCH --job-name=dm_register_analysis\n")
        out.write("#SBATCH --nodes=1\n")
        out.write("#SBATCH --ntasks=16\n")
        out.write("#SBATCH --mem-per-cpu=16g\n")
        out.write("#SBATCH --time=100:00:00\n")
        out.write("#SBATCH --mail-type=END\n")
        out.write("#SBATCH --no-requeue\n")

        cmd_upload = (
            f"module load java/11\n"
            f"export HPC_DM_UTILS={path_hpc_dme_clu}/utils\n"
            f"source $HPC_DM_UTILS/functions\n"
        )

        pi_name_path = f"{sample.get_attribute('PIfirstName')}_{sample.get_attribute('PIlastName')}"
        register_path = f"PI_Lab_{pi_name_path}"

        # Create PI Lab directory if it doesn't exist
        if not os.path.exists(register_path):
            os.mkdir(register_path)
            logging.info(f"Created PI Lab directory: {register_path}")
            create_pilab_json(register_path, sample.get_attribute('PIName'))
            cmd_upload += f"dm_register_collection {register_path}.metadata.json /FNL_SF_Archive/{register_path}\n"
        else:
            logging.info(f"PI Lab directory already exists: {register_path}")

        # Create Project directory
        register_path += f"/Project_{project_name}"
        if not os.path.exists(register_path):
            os.mkdir(register_path)
            logging.info(f"Created Project directory: {register_path}")
            create_project_json(
                register_path, 
                project_name, 
                sample.get_attribute('LabContact'), 
                sample.get_attribute('BAContact'), 
                sample.get_attribute('Project_date')
            )
            cmd_upload += f"dm_register_collection {register_path}.metadata.json /FNL_SF_Archive/{register_path}\n"
        else:
            logging.info(f"Project directory already exists: {register_path}")

        # Handle analysis files and metadata
        cmd_upload = handle_analysis_files(register_path, cmd_upload, sample.get_attribute('Application'))

        # Write the upload commands to the script
        out.write(cmd_upload)
        logging.info("Upload commands written to dm_register_analysis.sh")

def handle_analysis_files(register_path, cmd_upload, application):
    analysisfilelist = ''
    register_analysis_path = os.path.join(register_path, 'Analysis_combined')
    
    # Register Analysis_combined directory
    cmd_upload += f'dm_register_collection {register_analysis_path}.metadata.json /FNL_SF_Archive/{register_analysis_path}\n'
    
    if not os.path.exists(register_analysis_path):
        os.mkdir(register_analysis_path)
        logging.info(f"Created directory: {register_analysis_path}")
        today = datetime.today().strftime('%Y-%m-%d')
        create_flowcell_json(register_analysis_path, "Analysis", "nanopore", today, "combined", application)
    else:
        logging.info(f"Directory already exists: {register_analysis_path}")
        
    samplerun_in_metadata = read_sample_files()
    
    for samplerun in samplerun_in_metadata:
        samplename, runname = samplerun.split("|")
        register_sample_path = os.path.join(register_analysis_path, samplename)
        
        # Register each sample path
        cmd_upload += f'dm_register_collection {register_sample_path}.metadata.json /FNL_SF_Archive/{register_sample_path}\n'
        
        if not os.path.exists(register_sample_path):
            os.mkdir(register_sample_path)
            logging.info(f"Created directory: {register_sample_path}")
            create_sample_json(register_sample_path, samplename)
        else:
            logging.info(f"Directory already exists: {register_sample_path}")

        matrix_path = os.path.join(in_path, fproject, f"Sample_{samplename}", samplename)
        
        # Archive gene and transcript matrices
        folders_to_archive = [
            "gene_processed_feature_bc_matrix",
            "gene_raw_feature_bc_matrix",
            "transcript_processed_feature_bc_matrix",
            "transcript_raw_feature_bc_matrix"
        ]
        folder_paths = [os.path.join(matrix_path, folder) for folder in folders_to_archive]
        
        if all(os.path.exists(folder) for folder in folder_paths):
            tar_dest = os.path.join(register_sample_path, "gene_transcript_matrix.tar")
            try:
                with subprocess.Popen(
                    ["tar", "-cvhf", tar_dest, "-C", matrix_path] + folders_to_archive,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE
                ) as process:
                    stdout, stderr = process.communicate()
                    if process.returncode == 0:
                        logging.info(f"Successfully created tar file: {tar_dest}")
                    else:
                        logging.error(f"Error creating tar archive: {stderr.decode()}")
                        exit(1)
                analysisfilelist += f"{tar_dest}\n"
                create_object_json(register_sample_path, "gene_transcript_matrix.tar", samplename)
            except Exception as e:
                logging.error(f"Exception occurred while creating tar file for matrices: {e}")
        else:
            missing = [folder for folder in folder_paths if not os.path.exists(folder)]
            logging.warning(f"Missing folders for {samplename}: {', '.join(missing)}. Skipping archive creation.")
        
        # Create symlinks for BAM and BAI files
        bam_file = os.path.join(matrix_path, "tagged.bam")
        bai_file = os.path.join(matrix_path, "tagged.bam.bai")
        bam_symlink = os.path.join(register_sample_path, f"{samplename}.tagged.bam")
        bai_symlink = os.path.join(register_sample_path, f"{samplename}.tagged.bam.bai")
        
        if os.path.exists(bam_file) and os.path.exists(bai_file):
            try:
                os.symlink(bam_file, bam_symlink)
                os.symlink(bai_file, bai_symlink)
                logging.info(f"Created symlinks for BAM and BAI files in {register_sample_path}")
                analysisfilelist += f"{bam_symlink}\n{bai_symlink}\n"
                create_object_json(register_sample_path, f"{samplename}.tagged.bam", samplename)
                create_object_json(register_sample_path, f"{samplename}.tagged.bam.bai", samplename)
            except Exception as e:
                logging.error(f"Failed to create symlinks for BAM/BAI files: {e}")
        else:
            missing = [file for file in [bam_file, bai_file] if not os.path.exists(file)]
            logging.warning(f"Missing BAM or BAI files for {samplename}: {', '.join(missing)}. Skipping symlink creation.")

        # Tar and archive single-cell analysis outputs
        single_cell_files = [
            f'{samplename}_wf-single-cell-report.html', 'AllPCAPlot.pdf', 'TSNEPlots.pdf', 'UMAPPlots.pdf',
            'Filtered_Gene_Expression_Matrix.csv', 'Filtered_Normalized_Gene_Expression_Matrix.csv',
            'Gene_Expression_Matrix.csv', 'Number_of_cells_per_gene.csv', 'Number_of_genes_per_barcode.csv',
            'FilterNumbers.csv', 'annotations_BlueprintENCODE_general.csv', 'annotations_HPCA_general.csv',
            'ExpressionPlot.png', 'PreFilter_Gene_Plot.png', 'PostFilter_Gene_Plot.png',
            'GenesPerBarcodeEditedPlot.png', 'heatmap_blueprintencode.pdf', 'heatmap_hpca.pdf',
            'CellCycle.png', 'TSNEPlot_20PCs_0.1.png', 'TSNEPlot_20PCs_0.3.png', 'TSNEPlot_20PCs_0.6.png',
            'TSNEPlot_20PCs_0.8.png', 'top100markers_20_res0.1.csv', 'top100markers_20_res0.3.csv',
            'top100markers_20_res0.6.csv', 'top100markers_20_res0.8.csv', 'markers_res0.1.rds', 'markers_res0.3.rds',
            'markers_res0.6.rds', 'markers_res0.8.rds'
        ]
        sc_tar_dest = os.path.join(register_sample_path, "single_cell_analysis.tar")
        main_path = os.path.dirname(matrix_path)

        with tarfile.open(sc_tar_dest, "w") as tar:
            for file in single_cell_files:
                full_path = os.path.join(main_path, file)
                if os.path.exists(full_path):
                    tar.add(full_path, arcname=os.path.basename(file))
                else:
                    logging.warning(f"Single-cell file {file} not found. Skipping.")
        
        logging.info(f"Created single_cell_analysis.tar for {samplename}")
        analysisfilelist += f"{sc_tar_dest}\n"
        create_object_json(register_sample_path, "single_cell_analysis.tar", samplename)

        # Tar and archive IsoSeq analysis outputs
        isoseq_files = [
            f'{samplename}.annotated_transcripts.bed', f'{samplename}.sqanti_corrected.fasta',
            f'{samplename}.sqanti_classification_merge.txt', f'{samplename}.sqanti_SQANTI3_report.html',
            f'{samplename}.sqanti_RulesFilter_result_classification.txt',
            f'{samplename}.sqanti_filtering_reasons.txt', f'{samplename}.sqanti_junctions.txt',
            f'{samplename}.sqanti_SQANTI3_filter_report.pdf', f'{samplename}.sqanti_SQANTI3_report.html',
        ]
        isoseq_tar_dest = os.path.join(register_sample_path, "isoseq_analysis.tar")
        
        with tarfile.open(isoseq_tar_dest, "w") as tar:
            for file in isoseq_files:
                full_path = os.path.join(main_path, file)
                if os.path.exists(full_path):
                    tar.add(full_path, arcname=os.path.basename(file))
                else:
                    logging.warning(f"IsoSeq file {file} not found. Skipping.")
        
        logging.info(f"Created isoseq_analysis.tar for {samplename}")
        analysisfilelist += f"{isoseq_tar_dest}\n"
        create_object_json(register_sample_path, "isoseq_analysis.tar", samplename)

    # Save the analysis file list
    with open("FileList.txt", "w") as out_filelist:
        out_filelist.write(analysisfilelist)
    logging.info("Generated FileList.txt with paths of all archived files")
    
    return cmd_upload

def main():
    sample = load_metadata()
    config_working_directory(fproject, sample)

if __name__ == "__main__":
    main()
