# as reference, do not import directly, copy in the script and use as needed
# TODO: rewrite to shuttle about the cs instance as well to allow directly importing
# DONE: functions gauss, bimodal, get_threshold_of_scale_factor_bimodal, filter_particles_by_threshold, cs_import_particle_dataset can be used as imported functions

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.spatial.transform import Rotation as R
import ast

def create_next_workspace(cs, project, title):
    i=0
    j=0
    while j == 0:
        try:
            cs.find_workspace(project, f"W{i+1}")
            i+=1
        except:
            j = 1
    print(f"Next workspace in {project} is W{i+1}")
    ws = f"W{i+1}"
    ws = cs.create_workspace(project, title = title)
    return ws

def find_box_size(target_apix, tube_diameter = 376, factor = 1.6):
    good_box_sizes = [16, 20, 24, 28, 32, 36, 40, 48, 56, 60, 64, 72, 80, 84, 96, 100, 108, 112, 120, 128, 140, 144, 160, 168, 180, 192, 196, 200, 216, 224, 240, 252, 256, 280, 288, 300, 320, 324, 336, 360, 384, 392, 400, 420, 432, 448, 480, 500, 504, 512, 540, 560, 576, 588, 600, 640, 648, 672, 700, 720, 756, 768, 784, 800, 840, 864, 896, 900, 960, 972, 980, 1000, 1008, 1024, 1080, 1120, 1152, 1176, 1200, 1260, 1280, 1296, 1344, 1372, 1400, 1440, 1500, 1512, 1536, 1568, 1600, 1620, 1680, 1728, 1764, 1792, 1800, 1920, 1944, 1960, 2000]
    pixel_target = factor / target_apix * tube_diameter
    for box_size in good_box_sizes:
        if box_size < pixel_target:
            pass
        else:
            break
    return box_size


def log_append(logfile, message):
    with open(logfile, "a") as f:
        f.write(message + "\n")

def read_logfile(logfile_path):
    with open(logfile_path, "r") as f:
        for line in f:
            if line.startswith("Queued"):
                job, project, workspace, jobid, step = line.split()[1:6]
                yield job, project, workspace, jobid, step

def resume_pipeline(logfile_path, cs):
    pipeline = {}
    for jobNumber, project_id, workspace_id, job_id, step in read_logfile(logfile_path):
        #print(jobNumber, project_id, workspace_id, job_id, step)
        project = cs.find_project(project_id)
        workspace = project.find_workspace(workspace_id)
        job = project.find_job(job_id)
        print(f"Found job {jobNumber} in project {project.uid} workspace {workspace.uid} job {job.uid} step {step}.")
        pipeline[step] = job
        if job.status == "completed":
            print("Job is completed.")
        else:
            print(f"Job is {job.status}.")
            print(f"Continuing pipeline at stage {step}.")
            break
    return project, workspace, job, step, pipeline

def get_helical_refinement_results(job):
    a = job.model.dict()
    helical_order = a['spec']['outputs']['volume']['latest_summary_stats']['hsym_order']
    rise = a['spec']['outputs']['volume']['latest_summary_stats']['helical_rise_A']
    twist = a['spec']['outputs']['volume']['latest_summary_stats']['helical_twist_deg']
    return helical_order, rise, twist

def get_refinement_resolution(job):
    a = job.model.dict()
    resolution = a['spec']['outputs']['volume']['latest_summary_stats']['fsc_info_best']['radwn_final_A']

    return resolution

def get_particle_number(job):
    return len(job.load_output("particles"))


def get_refinement_nyquist(job):
    a = job.model.dict()
    psize = a['spec']['outputs']['particles']['summary']['blob/psize_A']
    nyquist = 2 * psize
    return nyquist

def set_selection_limits(job, dataset_fraction_per_class, nyquist_fraction):
    num_ptcls = get_particle_number(job)
    select_thresh_particles = int(dataset_fraction_per_class * num_ptcls)
    select_thresh_resolution = nyquist_fraction * get_refinement_nyquist(job)
    print(f"In job {job.uid} - {job.type}, selecting particles with a minimum of {select_thresh_particles} particles and a resolution of {select_thresh_resolution} A.")
    return select_thresh_particles, select_thresh_resolution

def cs_import_particle_dataset(project, workspace_uid, source_job, particles, title):
    """
    Args:
        project (_type_): _description_
        workspace_uid (_type_): _description_
        source_job (_type_): _description_
        particles (_type_): _description_
        title (_type_): _description_

    Returns:
        _type_: _description_
    """
    source_job.wait_for_done()
    for file in source_job.list_files():
        if file.endswith(".csg") and "particles" in file:
            
            output_name = file.split('.')[0] + "_uploaded"
            
            source_job.download_file(file, f"{file}")
            particles.save(f"{output_name}.cs")
            csg = []
            with open(file) as f:
                for line in f:
                    if "metafile:" not in line:
                        csg.append(line)
                    else:
                        csg.append(f"    metafile: \'>{output_name.split('.')[0]}.cs\'\n")
            csg = ''.join(csg)
            with open(f"{output_name}.csg", "w") as f:
                f.write(csg)
            
            project.upload(f"exports/groups/{source_job.uid}_particles/{output_name}.cs", f'{output_name}.cs', overwrite=True)
            project.upload(f"exports/groups/{source_job.uid}_particles/{output_name}.csg", f'{output_name}.csg', overwrite=True)
            
            upload = project.create_job(workspace_uid, "import_result_group",
                                    title=f"Import {title}",
                                    params={
                                        "blob_path": os.path.abspath(f"{project.dir()}/exports/groups/{source_job.uid}_particles/{output_name}.csg"),
                                    }
                                    )
            upload.queue()
            print(f"Importing {output_name}, {title} into project {project.uid} as {upload.uid} in workspace {workspace_uid}.")
            return upload
# filtered_particles = cs_import_particle_dataset(project, workspace_number, helix_refine, filtered_particles, "Filtered per-particle scale >= 0.8")



def gauss(x, mu, sigma, A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

def get_threshold_of_scale_factor_bimodal(job, sigma_mult=2, plot = False):
    """Deprecated. Use function from thresholding_functions.py instead.
    These accept cs job instance as input, leading to duplicated downloads from cs instance.
    New version shuttles about particles instead. 

    Args:
        job (_type_): _description_
        sigma_mult (int, optional): _description_. Defaults to 2.
        plot (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """
    particles = job.load_output("particles")
    df = pd.DataFrame(particles.rows())
    plt.hist(df['alignments3D/alpha'], bins = 100)
    y,x,_=plt.hist(df['alignments3D/alpha'], bins = 100)

    x=(x[1:]+x[:-1])/2 # for len(x)==len(y)

    #x, y inputs can be lists or 1D numpy arrays

    guess = (x[np.argmax(y)]/2, .2, np.max(y)/2, x[np.argmax(y)], .2, np.max(y))
    print(f"Two-gaussian fit guess: {guess}")

    params, cov = curve_fit(bimodal, x, y, guess)
    sigma=np.sqrt(np.diag(cov))
    x_fit = np.linspace(x.min(), x.max(), 500)
    #plot combined...
    plt.plot(x_fit, bimodal(x_fit, *params), color='red', lw=3, label='model')
    #...and individual Gauss curves
    plt.plot(x_fit, gauss(x_fit, *params[:3]), color='red', lw=1, ls="--", label='distribution 1')
    plt.plot(x_fit, gauss(x_fit, *params[3:]), color='red', lw=1, ls=":", label='distribution 2')
    #and the original data points if no histogram has been created before
    #plt.scatter(x, y, marker="X", color="black", label="original data")
    plt.xlabel("Scale factor")
    plt.ylabel("Frequency")
    plt.title(f"Filtering {job.project_uid} {job.uid} {job.type} by scale factor")
    
    if params[0] > params[3]:
        threshold = params[0] - 2 * params[1]
    else:
        threshold = params[3] - 2 * params[4]
    
    print(f"Threshold at {sigma_mult} sigma is {threshold}")
    
    plt.axvline(threshold, color='blue', lw=1, ls="--", label=f"{sigma_mult} sigma threshold")
    print(pd.DataFrame(data={'params': params, 'sigma': sigma}, index=bimodal.__code__.co_varnames[1:]))
    
    plt.legend()

    if plot:
        plt.show() 
        plt.close()
    else:
        plt.savefig(f"bimodal_fit_scale_{job.project_uid}_{job.uid}.png")
        plt.close()
    

    return threshold

def filter_particles_by_threshold(job, threshold, field, plot=False):
    """Deprecated. Use function from thresholding_functions.py instead.
    These accept cs job instance as input, leading to duplicated downloads from cs instance.
    New version shuttles about particles instead. 

    Args:
        job (_type_): _description_
        sigma_mult (int, optional): _description_. Defaults to 2.
        plot (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
        
    For scale, field = 'alignments3D/alpha' """
    
    particles = job.load_output("particles")
    filtered_particles = particles.query(
    lambda x: x[field] > threshold
    )
    refused_particles = particles.query(
        lambda x: x[field] < threshold
    )
    """if plot:
        df = pd.DataFrame(filtered_particles.rows())
        plt.hist(df[field], bins = 100)
        plt.xlabel(field)
        plt.ylabel("Frequency")
        plt.title(f"Filtered particles by {field} > {threshold}")
        plt.annotate(f"Particles remaining: {len(filtered_particles)}", (0.6, 0.9), xycoords='axes fraction')
        plt.savefig(f"filtered_particles_{job.uid}.png")"""
    return filtered_particles, refused_particles

def filter_particles_by_tilt_threshold(job, threshold1, threshold2):
    """ For thresholding on pose.
    Deprecated. Use function from thresholding_functions.py instead.
    These accept cs job instance as input, leading to duplicated downloads from cs instance.
    New version shuttles about particles instead. 
    """
    field = "tilt"
    t1 = threshold1
    t2 = threshold2
    particles = job.load_output("particles")
    filtered_particles = particles.query(
    lambda x: R.from_rotvec(x["alignments3D/pose"]).as_euler("ZYZ", degrees=True)[1] > t1 and R.from_rotvec(x["alignments3D/pose"]).as_euler("ZYZ", degrees=True)[1] < t2
    )
    refused_particles = particles.query(
        lambda x: R.from_rotvec(x["alignments3D/pose"]).as_euler("ZYZ", degrees=True)[1] < t1 or R.from_rotvec(x["alignments3D/pose"]).as_euler("ZYZ", degrees=True)[1] > t2
    )
    return filtered_particles, refused_particles

def get_threshold_of_tilt_bimodal(job, sigma_mult = 2, plot=False):
    """Deprecated. Use function from thresholding_functions.py instead.
    These accept cs job instance as input, leading to duplicated downloads from cs instance.
    New version shuttles about particles instead. 
    """
    particles = job.load_output("particles")
    eulers = R.from_rotvec(particles["alignments3D/pose"]).as_euler("ZYZ", degrees=True)
    #psi = eulers[:, 2]
    tilt = eulers[:, 1] 
    #phi = eulers[:, 0]
    y,x,_=plt.hist(tilt, bins = 100)
    x=(x[1:]+x[:-1])/2 # for len(x)==len(y)

    #x, y inputs can be lists or 1D numpy arrays

    guess = (x[np.argmax(y)], 25, np.max(y)/10, x[np.argmax(y)], 5, np.max(y))
    print(f"Two-gaussian fit guess: {guess}")
    try:
        print("Calculating bimodal fit")
        params, cov = curve_fit(bimodal, x, y, guess)
        sigma=np.sqrt(np.diag(cov))
        x_fit = np.linspace(x.min(), x.max(), 500)
        #plot combined...
        plt.plot(x_fit, bimodal(x_fit, *params), color='red', lw=3, label='model')
        #...and individual Gauss curves
        plt.plot(x_fit, gauss(x_fit, *params[:3]), color='red', lw=1, ls="--", label='distribution 1')
        plt.plot(x_fit, gauss(x_fit, *params[3:]), color='red', lw=1, ls=":", label='distribution 2')
        #and the original data points if no histogram has been created before
        #plt.scatter(x, y, marker="X", color="black", label="original data")
        plt.xlabel("Tilt")
        plt.ylabel("Frequency")
        #plt.title(f"Filtering {job.uid} {job.type} by scale factor")

        if params[0] > params[3]:
            threshold1 = params[0] - 2 * params[1]
            threshold2 = params[0] + 2 * params[1]
        else:
            threshold1 = params[3] - 2 * params[4]
            threshold2 = params[3] + 2 * params[4]
    except:
        print("Bimodal fit failed. Trying with single gaussian.")
        guess = guess[3:]
        params, cov = curve_fit(gauss, x, y, guess)
        sigma=np.sqrt(np.diag(cov))
        x_fit = np.linspace(x.min(), x.max(), 500)
        #plot combined...
        #...and individual Gauss curves
        plt.plot(x_fit, gauss(x_fit, *params[:3]), color='red', lw=1, ls="--", label='distribution 1')
        #and the original data points if no histogram has been created before
        #plt.scatter(x, y, marker="X", color="black", label="original data")
        plt.xlabel("Tilt")
        plt.ylabel("Frequency")
        #plt.title(f"Filtering {job.uid} {job.type} by scale factor")

        threshold1 = params[0] - 2 * params[1]
        threshold2 = params[0] + 2 * params[1]
    

    print(f"Threshold at {sigma_mult} sigma is {threshold1} - {threshold2}")

    plt.axvline(threshold1, color='blue', lw=1, ls="--", label=f"{sigma_mult} sigma threshold")
    plt.axvline(threshold2, color='blue', lw=1, ls="--")
    #print(pd.DataFrame(data={'params': params, 'sigma': sigma}, index=bimodal.__code__.co_varnames[1:]))

    plt.legend()
    if plot:
        plt.show() 
        plt.close()
    else:
        plt.savefig(f"bimodal_fit_tilt_{job.project_uid}_{job.uid}.png")
        plt.close()
    return threshold1, threshold2   