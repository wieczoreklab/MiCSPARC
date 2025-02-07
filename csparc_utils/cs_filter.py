

from cryosparc.tools import CryoSPARC
from thresholding_functions import gauss, bimodal, get_threshold_of_scale_factor_bimodal, get_threshold_of_tilt_bimodal, filter_particles_by_threshold, filter_particles_by_tilt_threshold, cs_import_particle_dataset
import os
import yaml
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--project', required=True, help='Project name')
    parser.add_argument('--workspace', required=True, help='Workspace name')
    parser.add_argument('--job', required=True, help='Job name')
    parser.add_argument('--filter', required=True, help='Type of filter. Options: scale, tilt, scale_tilt, scale_plain')
    parser.add_argument('--threshold', required=False, help='Threshold value for scale_plain filter')
    parser.add_argument('--threshold_sigma', required=False, default=2, help='Sigma for threshold for scale and tilt filters, default 2')
    parser.add_argument('--keep_refused', required=False, default=False, help='Keep refused particles from the last filter. Useful for scale filter.')
    args = parser.parse_args()

    homedir = os.getenv("HOME")
    with open(f"{homedir}/cs_config.yml", "r") as f:
        cs_conf = yaml.safe_load(f.read())


    cs = CryoSPARC(
        license = cs_conf["license"],
        host = cs_conf['host'],
        base_port = cs_conf['base_port'],
        email = cs_conf['email'],
        password = cs_conf['password']
    )
    assert cs.test_connection()
    

    
    project = args.project
    workspace = args.workspace
    job = args.job
    threshold_sigma = float(args.threshold_sigma)
    filter = args.filter
    keep_refused = args.keep_refused
    project = cs.find_project(project)
    job = project.find_job(job)
    particles = job.load_output("particles")
    
    
    if filter == 'scale':
        threshold = get_threshold_of_scale_factor_bimodal(job, particles, sigma_mult=threshold_sigma)
        particles, refused = filter_particles_by_threshold(particles, threshold, "alignments3D/alpha")
        uploaded = cs_import_particle_dataset(project, 
                                              workspace, 
                                              job, 
                                              particles, 
                                              f"scale over {threshold:2f} filtered particles from {job.uid}")
    elif filter == 'scale_plain':
        threshold = float(args.threshold)  
        particles, refused = filter_particles_by_threshold(particles, threshold, "alignments3D/alpha")
        uploaded = cs_import_particle_dataset(project, 
                                              workspace, 
                                              job, 
                                              particles, 
                                              f"scale over {threshold:2f} filtered particles from {job.uid}")
    elif filter == 'tilt':
        t1, t2 = get_threshold_of_tilt_bimodal(job, particles, sigma_mult=threshold_sigma)
        particles, refused = filter_particles_by_tilt_threshold(particles, t1, t2)
        uploaded = cs_import_particle_dataset(project, 
                                              workspace, 
                                              job, 
                                              particles, 
                                              f"tilt between {t1:2f} and {t2:2f} filtered particles from {job.uid}")
        
    elif filter == 'scale_tilt':
        threshold = get_threshold_of_scale_factor_bimodal(job, particles, sigma_mult=threshold_sigma)
        particles, refused = filter_particles_by_threshold(particles, threshold, "alignments3D/alpha")
        uploaded = cs_import_particle_dataset(project, 
                                              workspace, 
                                              job, 
                                              particles, 
                                              f"scale over {threshold:2f} filtered particles from {job.uid}")
        t1, t2 = get_threshold_of_tilt_bimodal(job, particles, sigma_mult=threshold_sigma)
        particles, refused = filter_particles_by_tilt_threshold(particles, t1, t2)
        uploaded = cs_import_particle_dataset(project, 
                                              workspace, 
                                              uploaded, 
                                              particles, 
                                              f"tilt between {t1:2f} and {t2:2f} filtered particles from {uploaded.uid}")
    else:
        raise ValueError('Unknown threshold type')
    
    print(f"{len(particles)} particles remaining.")
    print(f"{len(refused)} particles refused.")
    print(f"Filtered particles uploaded to {uploaded.uid}")

    if keep_refused:
        print(f"Uploading refused particles")
        upload_refuse = cs_import_particle_dataset(project, 
                                              workspace, 
                                              job, 
                                              refused, 
                                              f"Refused particles from {job.uid}")
        print(f"Refused particles uploaded to {upload_refuse.uid}")
