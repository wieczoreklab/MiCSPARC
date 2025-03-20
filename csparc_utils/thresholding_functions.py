# filtering functions rewrite to shuffle about particles rather than job instances
# for speed with large datasets

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.spatial.transform import Rotation as R
import ast

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

def get_threshold_of_scale_factor_bimodal(source_job, particles, sigma_mult=2, plot = False):
    project = source_job.project_uid
    job = source_job.uid
    jtype = source_job.type
    
    
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
    plt.title(f"Filtering {project} {job} {jtype} by scale factor")
    
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
        plt.savefig(f"bimodal_fit_scale_{project}_{job}.png")
        plt.close()
    return threshold

def get_threshold_of_tilt_bimodal(source_job, particles, sigma_mult = 2, plot=False):
    project = source_job.project_uid
    job = source_job.uid
    jtype = source_job.type
    
    eulers = R.from_rotvec(particles["alignments3D/pose"]).as_euler("ZYZ", degrees=True)
    #psi = eulers[:, 2]
    tilt = eulers[:, 1] 
    #phi = eulers[:, 0]
    y,x,_=plt.hist(tilt, bins = 100)
    x=(x[1:]+x[:-1])/2 # for len(x)==len(y) #x, y inputs can be lists or 1D numpy arrays
    guess = (x[np.argmax(y)], 25, np.max(y)/10, x[np.argmax(y)], 5, np.max(y)) # One wide one narrow gaussian with shared centre
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
        plt.plot(x_fit, gauss(x_fit, *params[:3]), color='red', lw=1, ls="--", label='model')
        #and the original data points if no histogram has been created before
        #plt.scatter(x, y, marker="X", color="black", label="original data")
        plt.xlabel("Tilt")
        plt.ylabel("Frequency")
        #plt.title(f"Filtering {job.uid} {job.type} by scale factor")

        threshold1 = params[0] - 2 * params[1]
        threshold2 = params[0] + 2 * params[1]
    

    print(f"Threshold at {sigma_mult} sigma is {threshold1} - {threshold2}")
    
    plt.title(f"Filtering {project} {job} {jtype} by tilt")

    plt.axvline(threshold1, color='blue', lw=1, ls="--", label=f"{sigma_mult} sigma threshold")
    plt.axvline(threshold2, color='blue', lw=1, ls="--")
    #print(pd.DataFrame(data={'params': params, 'sigma': sigma}, index=bimodal.__code__.co_varnames[1:]))

    plt.legend()
    if plot:
        plt.show() 
        plt.close()
    else:
        plt.savefig(f"bimodal_fit_tilt_{project}_{job}.png")
        plt.close()
        
    return threshold1, threshold2   

def filter_particles_by_threshold(particles, threshold, field, plot=False):
    """ For scale, field = 'alignments3D/alpha' """
    filtered_particles = particles.query(
    lambda x: x[field] > threshold
    )
    refused_particles = particles.query(
        lambda x: x[field] < threshold
    )
    return filtered_particles, refused_particles

def filter_particles_by_tilt_threshold(particles, threshold1, threshold2):
    """ For thresholding on pose """
    field = "tilt"
    t1 = threshold1
    t2 = threshold2
    filtered_particles = particles.query(
    lambda x: R.from_rotvec(x["alignments3D/pose"]).as_euler("ZYZ", degrees=True)[1] > t1 and R.from_rotvec(x["alignments3D/pose"]).as_euler("ZYZ", degrees=True)[1] < t2
    )
    refused_particles = particles.query(
        lambda x: R.from_rotvec(x["alignments3D/pose"]).as_euler("ZYZ", degrees=True)[1] < t1 or R.from_rotvec(x["alignments3D/pose"]).as_euler("ZYZ", degrees=True)[1] > t2
    )
    return filtered_particles, refused_particles

