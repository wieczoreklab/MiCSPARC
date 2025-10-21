#!/data/user/filipc_p/miniforge3/envs/micsparc/bin/python
import os
import multiprocessing as mp
import click
import numpy as np
from tqdm import tqdm
import pandas as pd
import yaml
import time
import ast
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import differential_evolution
from scipy.spatial.transform import Rotation as R
from scipy import stats
from cryosparc.tools import CryoSPARC
from cryosparc.dataset import Dataset
from csparc_unify_phi import phi_process_mic
from cryosparc_functions import create_next_workspace, find_box_size, log_append, read_logfile, resume_pipeline, get_helical_refinement_results, get_refinement_resolution, get_particle_number, get_refinement_nyquist, set_selection_limits
from csparc_extrapolate_filaments_hmh import extrapolate_filaments
from csparc_assign_pfns import class_posterior


if __name__ == "__main__":
    #Instance configuration
    
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


    logfile_path = f"autopick.log"

    #Lane configuration

    cpu_lane = cs_conf.get('cpu_lane', "cpu-daily")
    gpu_lane = cs_conf.get('gpu_lane', "gpu")
    cpu_highmem_lane = cs_conf.get('cpu_highmem_lane', "cpu-highmem")
    gpu_highmem_lane = cs_conf.get('gpu_highmem_lane', "gpu-highmem")
    gpu_highmem_short_lane = cs_conf.get('gpu_highmem_short_lane', "gpu-highmem-short")
    gpu_veryhighmem_lane = cs_conf.get('gpu_veryhighmem_lane', "gpu-node")
    gpu_veryhighmem_short_lane = cs_conf.get('gpu_veryhighmem_short_lane', "gpu-node-short")
    models = cs_conf.get('models', "/data/project/bio/steinmetz/software/microtubule-processing/models/kinesin/*aligned.mrc")

    scratch = cs_conf.get('scratch', True)

    if os.path.isfile(logfile_path):
        project, ws, job, step, pipeline = resume_pipeline(logfile_path, cs)
        with open(logfile_path, "r") as f:
            for line in f:
                if line.startswith(" Segment length"):
                    segment_length = int(line.split(":")[1])
                if line.startswith(" Tube diameter"):
                    tube_diameter = int(line.split(":")[1])
                if line.startswith(" Particle padding factor"):
                    particle_padding_factor = float(line.split(":")[1])
                if line.startswith(" Apix"):
                    apix_original = float(line.split(":")[1])
                if line.startswith("Valid classes"):
                    valid_classes = line.split(":")[1]
                    valid_classes = ast.literal_eval(valid_classes)
                    print(f"Valid classes: {valid_classes}")
                    
        print(f"Resuming pipeline at step {step}.")
        print(f"Segment length: {segment_length}")
        print(f"Tube diameter: {tube_diameter}")
        print(f"Particle padding factor: {particle_padding_factor}")
        print(f"Apix: {apix_original}")
        print(f"Project: {project.uid}")
        print(f"Workspace: {ws.uid}")
        resume = True
        
    else:
        proj = input("Which project do you want to use? [P1, P2, etc.]:")
        import_jobid = input("What is the ID of your import job [J1, J2, etc.]:")
        project = cs.find_project(proj)
        importjob = project.find_job(import_jobid)
        project_directory = project.dir()
        os.chdir(project_directory)
        
        #Autopicking settings

        if os.path.isfile(f"cs_autopick.yml"):
            with open(f"cs_autopick.yml", "r") as f:
                pick_settings = yaml.safe_load(f.read())
                segment_length = pick_settings.get('segment_length', 82)
                tube_diameter = pick_settings.get('tube_diameter', 376)
                particle_padding_factor = pick_settings.get('particle_padding_factor', 1.7)
                models = pick_settings.get('models', models)
        else:
            segment_length = int(input("What is the segment length in Ångstrom? [default: 82]:") or 82)
            tube_diameter = int(input("What is the tube diameter in Ångstrom? [default: 376]:") or 376)
            particle_padding_factor = float(input("What is the particle padding factor? [default: 1.7]:") or 1.7)
            with open(f"cs_autopick.yml", "w") as f:
                yaml.dump({"segment_length": segment_length, "tube_diameter": tube_diameter, "particle_padding_factor": particle_padding_factor}, f)

        pipeline = {}
        pipeline["import"] = importjob
        importjob.wait_for_done()
        ij = importjob.load_output("imported_movies")
        apix_original = ij['movie_blob/psize_A'][0]
        pixel_target_original = particle_padding_factor / apix_original * tube_diameter

        ws = create_next_workspace(cs, proj, "Autopicking")
        
        log_append(logfile_path, f"Starting autopicking for project {proj} with import job {import_jobid}")
        autopick_settings = f"Picking will be done with the following settings: \n Segment length: {segment_length} \n Tube diameter: {tube_diameter} \n Particle padding factor: {particle_padding_factor} \n Apix: {apix_original} \n"
        log_append(logfile_path, autopick_settings)
        resume = False
        
    microtubule_params_theory = {
    '11-3': {
        'rise': 11.2,
        'twist': -32.7
    },
    '12-3': {
        'rise': 10.3,
        'twist': -30.0
    },
    '13-3': {
        'rise': 9.46,
        'twist': -27.7
    },
    '14-3': {
        'rise': 8.79,
        'twist': -25.7
    },
    '15-4': {
        'rise': 10.9,
        'twist': -24
    },
    '16-4': {
        'rise': 10.3,
        'twist': -22.5
    }
    }
    
    jobcounter = 1
    if resume == False or step == "motioncor_small":
        resume = False
        
        pipeline['motioncor_small'] = project.create_job(ws.uid, 
                                    "patch_motion_correction_multi", 
                                    connections={
                                        "movies": (importjob.uid, "imported_movies")
                                    }, 
                                    params = {
                                        "random_num" : 500,
                                        "output_f16" : True
                                    })
        pipeline['motioncor_small'].queue(lane = gpu_highmem_short_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['motioncor_small'].uid} motioncor_small")
        jobcounter += 1
        
    if resume == False or step == "motioncor_full":
        resume = False
        
        pipeline['motioncor_full'] = project.create_job(ws.uid, 
                                    "patch_motion_correction_multi", 
                                    connections={
                                        "movies": (importjob.uid, "imported_movies")
                                    }, 
                                    params = {
                                        "output_f16": True
                                    })
        pipeline['motioncor_full'].queue(lane = gpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['motioncor_full'].uid} motioncor_full")
        jobcounter += 1
        print("Motion correction submitted.")

    if resume == False or step == "ctffind_small":
        resume = False
        
        pipeline['ctffind_small'] = project.create_job(ws.uid, 
                                    "patch_ctf_estimation_multi", 
                                    connections={
                                        "exposures": (pipeline['motioncor_small'].uid, "micrographs")
                                    }, 
                                    params = {
                                        
                                    })
        pipeline['ctffind_small'].queue(lane = gpu_highmem_short_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['ctffind_small'].uid} ctffind_small")
        jobcounter += 1
        

    if resume == False or step == "ctffind_full":
        resume = False
        
        pipeline['ctffind_full'] = project.create_job(ws.uid,
                                            "patch_ctf_estimation_multi",
                                            connections={
                                                "exposures": (pipeline['motioncor_full'].uid, "micrographs")
                                            },
                                            params = {
                                            }
                                            )
        pipeline['ctffind_full'].queue(lane = gpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['ctffind_full'].uid} ctffind_full")
        jobcounter += 1
        print("CTF estimation submitted.")

    if resume == False or step == "curate_full":
        resume = False
        
        pipeline['curate_full'] = project.create_job(ws.uid, 
                                "curate_exposures_v2", 
                                title="Filter to 2-4.5A CTF fits",
                                connections={
                                    "exposures": (pipeline['ctffind_full'].uid, "exposures")
                                }, 
                                params = {
                                    "auto_threshold_ctf_fit_to_A": "2, 4.5"
                                })
        pipeline['curate_full'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['curate_full'].uid} curate_full")
        jobcounter += 1
        
    if resume == False or step == "tracer1":
        resume = False
        
        pipeline['tracer1'] = project.create_job(ws.uid,
                                "filament_tracer_gpu",
                                title = "First filament trace",
                                params={
                                    "filament_diameter": tube_diameter,
                                    "min_distance": segment_length/tube_diameter,
                                    "diameter": tube_diameter - 76,
                                    "diameter_max": tube_diameter + 76,
                                    "num_process": 500
                                }
    )
        pipeline['tracer1'].connect("micrographs", pipeline['ctffind_small'].uid, "exposures")
        pipeline['tracer1'].queue(lane= gpu_highmem_short_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['tracer1'].uid} tracer1")
        jobcounter += 1

    stage_apix = 4    

    if resume == False or step == "extract1":
        resume = False
        
        print(f"Stage 1, processing at {stage_apix} A/px")
        pipeline['extract1'] = project.create_job(ws.uid, 
                                    "extract_micrographs_cpu_parallel", 
                                    title="1st Extract and bin to 4 A / px",
                                    connections={
                                        "micrographs": (pipeline['tracer1'].uid, "micrographs"),
                                        "particles": (pipeline['tracer1'].uid, "particles")
                                    }, 
                                    params = {
                                        "box_size_pix" : find_box_size(apix_original, tube_diameter=tube_diameter, factor=particle_padding_factor),
                                        "bin_size_pix" : find_box_size(stage_apix, tube_diameter=tube_diameter, factor=particle_padding_factor),
                                    })
        pipeline['extract1'].queue(lane = cpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['extract1'].uid} extract1")
        jobcounter += 1

    if resume == False or step == "classify_blob_picks":
        resume = False
        
        pipeline['classify_blob_picks'] = project.create_job(ws.uid,
                                                    "class_2D_new",
                                                    title="1st picks classify",
                                                    connections={"particles": (pipeline['extract1'].uid, "particles")},
                                                    params={
                                                        "class2D_K": 10,
                                                        "class2D_estimate_in_plane_pose": True,
                                                        "compute_use_ssd": scratch,
                                                        },
        )
        pipeline['classify_blob_picks'].queue(lane= gpu_highmem_short_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['classify_blob_picks'].uid} classify_blob_picks")
        jobcounter += 1

    if resume == False or step == "select2d_1":
        resume = False
        
        pipeline["classify_blob_picks"].wait_for_done()
        
        particle_number = get_particle_number(pipeline['classify_blob_picks'])
        if particle_number < 10000:
            nq_frac = 5
        else:
            nq_frac = 2.5
        
        thresh_particles, thresh_resolution = set_selection_limits(
            pipeline['classify_blob_picks'], 
            dataset_fraction_per_class = 0.1, 
            nyquist_fraction = nq_frac)
        
        pipeline['select2d_1'] = project.create_job(ws.uid,
                                        "select_2D",
                                        title = f"Select more than {thresh_particles} ptcl, better than {thresh_resolution} A resolution",
                                        params={
                                            "particle_count_above": thresh_particles,
                                            "resolution_better_than": thresh_resolution
                                            },
                                        connections={
                                            "particles": (pipeline['classify_blob_picks'].uid, "particles"),
                                            "templates": (pipeline['classify_blob_picks'].uid, "class_averages"),
                                            })
        pipeline['select2d_1'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['select2d_1'].uid} select2d_1")
        jobcounter += 1
        
    if resume == False or step == "template_tracer1":
        resume = False
        
        pipeline['template_tracer1'] = project.create_job(
                                            ws.uid, 
                                            "filament_tracer_gpu",
                                            title = "1st template tracer",
                                            params={"filament_diameter": tube_diameter,
                                                    "min_distance": segment_length/tube_diameter,
                                                    "min_filament_length": 1,
                                                    "branch_radius_remove": 0.5,
                                                    "max_prune_dist": 0.5,
                                                    "num_process": 500
                                                    }
                                        )
        pipeline['template_tracer1'].connect("micrographs", pipeline['ctffind_small'].uid, "exposures")
        pipeline['template_tracer1'].connect("templates", pipeline['select2d_1'].uid, "templates_selected")
        pipeline['template_tracer1'].queue(lane= gpu_highmem_short_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['template_tracer1'].uid} template_tracer1")
        jobcounter += 1

    if resume == False or step == "extract2":
        resume = False
        pipeline['extract2'] = project.create_job(ws.uid,
                                        "extract_micrographs_cpu_parallel",
                                        title="2nd Extract and bin to 4 A/px",
                                        connections={
                                            "micrographs": (pipeline['template_tracer1'].uid, "micrographs"),
                                            "particles": (pipeline['template_tracer1'].uid, "particles")
                                        },
                                        params={
                                            "box_size_pix": find_box_size(apix_original, tube_diameter=tube_diameter, factor=particle_padding_factor),
                                            "bin_size_pix": find_box_size(stage_apix, tube_diameter=tube_diameter, factor=particle_padding_factor),
                                        })
        pipeline['extract2'].queue(lane= cpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['extract2'].uid} extract2")
        jobcounter += 1
        
    if resume == False or step == "classify_first_picks":
        resume = False
        
        pipeline['classify_first_picks'] = project.create_job(ws.uid,
                                    "class_2D_new",
                                    title="Classify first template picks",
                                    connections={"particles": (pipeline['extract2'].uid, "particles")},
                                    params={"class2D_K": 50,
                                            "class2D_estimate_in_plane_pose": True,
                                            "compute_use_ssd": scratch,
                                        },
                                )
        pipeline['classify_first_picks'].queue(lane=gpu_highmem_short_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['classify_first_picks'].uid} classify_first_picks")
        jobcounter += 1
        
    if resume == False or step == "select2d_2":
        resume = False
        
        pipeline['classify_first_picks'].wait_for_done()
        
        particle_number = get_particle_number(pipeline['classify_first_picks'])
        
        if particle_number < 10000:
            nq_frac = 5
        else:
            nq_frac = 2.2
        
        thresh_particles, thresh_resolution = set_selection_limits(
            pipeline['classify_first_picks'], 
            dataset_fraction_per_class = 0.05, 
            nyquist_fraction = nq_frac)
        
        pipeline['select2d_2'] = project.create_job(ws.uid, 
                                                    "select_2D",
                                                    title=f"Select more than {thresh_particles} ptcl, better than {thresh_resolution} A resolution",
                                    connections={
            "particles": (pipeline['classify_first_picks'].uid, "particles"),
            "templates": (pipeline['classify_first_picks'].uid, "class_averages"),
        })
        pipeline['select2d_2'].queue()
        pipeline['select2d_2'].wait_for_status("waiting")
        class_info = pipeline['select2d_2'].interact("get_class_info")
        for c in class_info:
            if 1.0 < c["res_A"] < thresh_resolution and c["num_particles_total"] > thresh_particles:
                pipeline['select2d_2'].interact(
                    "set_class_selected",
                    {
                        "class_idx": c["class_idx"],
                        "selected": True,
                    },
                )
        pipeline['select2d_2'].interact("finish")
        pipeline['select2d_2'].wait_for_done()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['select2d_2'].uid} select2d_2")
        jobcounter += 1

    if resume == False or step == "template_tracer2":
        resume = False
        
        pipeline['template_tracer2'] = project.create_job(
            ws.uid,
            "filament_tracer_gpu",
            title = "Second template trace, full dataset",
            params={"filament_diameter": tube_diameter,
                    "min_distance": segment_length/tube_diameter,
                    "min_filament_length": 1,
                    "branch_radius_remove": 0.5,
                    "max_prune_dist": 0.5
                    }
        )
        pipeline['template_tracer2'].connect("micrographs", pipeline['curate_full'].uid, "exposures_accepted")
        pipeline['template_tracer2'].connect("templates", pipeline['select2d_2'].uid, "templates_selected")
        pipeline['template_tracer2'].queue(lane=gpu_lane)
        
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['template_tracer2'].uid} template_tracer2")
        jobcounter += 1
        
    if resume == False or step == "inspect1":
        resume = False
        
        pipeline['inspect1'] = project.create_job(
            ws.uid,
            "inspect_picks_v2",
            connections = {
                "micrographs" : (pipeline['curate_full'].uid, "exposures_accepted"),
                "particles" : (pipeline['template_tracer2'].uid, "particles")
            }
        )
        pipeline['inspect1'].queue()

        pipeline['template_tracer2'].wait_for_done()
            
        picks = pipeline['template_tracer2'].load_output("particles")

        picks_pd = pd.DataFrame(picks.rows())

        curvature_mean = np.mean(picks_pd['filament/curvature'])
        curvature_sd = np.std(picks_pd['filament/curvature'])
        sinuosity_mean = np.mean(picks_pd['filament/sinuosity'])
        sinuosity_sd = np.std(picks_pd['filament/sinuosity'])
        pick_power_mean = np.mean(picks_pd['pick_stats/power'])
        pick_power_sd = np.std(picks_pd['pick_stats/power'])
        ncc_score_mean = np.mean(picks_pd['pick_stats/ncc_score'])
        ncc_score_sd = np.std(picks_pd['pick_stats/ncc_score'])
        print(f"Second filament trace results. Filament curvature is {curvature_mean} +- {curvature_sd}.")
        print(f"Sinuosity is {sinuosity_mean} +- {sinuosity_sd}.")
        print(f"Pick power is {pick_power_mean} +- {pick_power_sd}.")
        print(f"NCC Score is {ncc_score_mean} +- {ncc_score_sd}.")

        thresh_curv = curvature_mean + curvature_sd
        thresh_sin = sinuosity_mean + sinuosity_sd
        thresh_pick = pick_power_mean - pick_power_sd
        thresh_pick_max  = pick_power_mean +2* pick_power_sd
        thresh_ncc = ncc_score_mean - ncc_score_sd

        pipeline['inspect1'].wait_for_status("waiting")
        pipeline['inspect1'].interact(
                        "set_thresholds",
                        {'curv_thresh': float(thresh_curv),
                        'sinu_thresh': float(thresh_sin),
                        'lpower_thresh_min': int(thresh_pick),
                        'lpower_thresh_max': int(thresh_pick_max),
                        'ncc_score_thresh': float(thresh_ncc)},
                        )
        pipeline['inspect1'].interact("shutdown_interactive")
        pipeline['inspect1'].wait_for_done()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['inspect1'].uid} inspect1")
        jobcounter += 1
        
    if resume == False or step == "extract3":
        resume = False
        
        pipeline['extract3'] = project.create_job(ws.uid,
                                        "extract_micrographs_cpu_parallel",
                                        title=f"3rd Extract and bin to {stage_apix} A/px",
                                        connections={
                                            "micrographs": (pipeline['inspect1'].uid, "micrographs"),
                                            "particles": (pipeline['inspect1'].uid, "particles")
                                        },
                                        params={
                                            "box_size_pix": find_box_size(apix_original, tube_diameter=tube_diameter, factor=particle_padding_factor),
                                            "bin_size_pix": find_box_size(stage_apix, tube_diameter=tube_diameter, factor=particle_padding_factor),
                                        })
        pipeline['extract3'].queue(lane=cpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['extract3'].uid} extract3")
        jobcounter += 1
        
    if resume == False or step == "classify_good_template_picks":
        resume = False
        
        pipeline['classify_good_template_picks'] = project.create_job(ws.uid,
                                                        "class_2D_new",
                                                            title="Good template picks classify 1",
                                                            connections={"particles": (pipeline['extract3'].uid, "particles")},
                                                            params={
                                                                "class2D_K": 100,
                                                                "class2D_estimate_in_plane_pose": True,
                                                                "class2D_num_full_iter": 2,
                                                                "compute_use_ssd": scratch,
                                                            },
        )
        pipeline['classify_good_template_picks'].queue(lane=gpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['classify_good_template_picks'].uid} classify_good_template_picks")
        jobcounter += 1
        
    if resume == False or step == "select2d_4":
        resume = False
        
        pipeline['classify_good_template_picks'].wait_for_done()
        
        particle_number = get_particle_number(pipeline['classify_good_template_picks'])
        if particle_number < 10000:
            nq_frac = 3
        else:
            nq_frac = 2.1
        
        thresh_particles, thresh_resolution = set_selection_limits(
            pipeline['classify_good_template_picks'], 
            dataset_fraction_per_class = 0.03, 
            nyquist_fraction = nq_frac)
        
        pipeline['select2d_4'] = project.create_job(ws.uid, 
                                        "select_2D",
                                        title = f"Select more than {thresh_particles} ptcl, better than {thresh_resolution} A resolution",
                                        params={
                                            "particle_count_above": thresh_particles,
                                            "resolution_better_than": thresh_resolution
                                            },
                                        connections={
                                            "particles": (pipeline['classify_good_template_picks'].uid, "particles"),
                                            "templates": (pipeline['classify_good_template_picks'].uid, "class_averages"),
                                            })
        pipeline['select2d_4'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['select2d_4'].uid} select2d_4")
        jobcounter += 1
        
    if resume == False or step == "classify_good_template_picks2":
        resume = False
        
        pipeline['classify_good_template_picks2'] = project.create_job(ws.uid,
                                                        "class_2D_new",
                                                            title="Good template picks classify 2",
                                                            connections={"particles": (pipeline['select2d_4'].uid, "particles_selected")},
                                                            params={
                                                                "class2D_K": 50,
                                                                "class2D_estimate_in_plane_pose": True,
                                                                "class2D_num_full_iter": 2,
                                                                "class2D_sigma_init_iter": 200,
                                                                "compute_use_ssd": scratch,
                                                            },
        )
        pipeline['classify_good_template_picks2'].queue(lane=gpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['classify_good_template_picks2'].uid} classify_good_template_picks2")
        jobcounter += 1

    if resume == False or step == "select2d_5":
        resume = False
        
        pipeline['classify_good_template_picks2'].wait_for_done()
        
        particle_number = get_particle_number(pipeline['classify_good_template_picks2'])
        if particle_number < 10000:
            nq_frac = 3
        else:
            nq_frac = 2.1
        
        thresh_particles, thresh_resolution = set_selection_limits(
            pipeline['classify_good_template_picks2'], 
            dataset_fraction_per_class = 0.02, 
            nyquist_fraction = nq_frac)
        
        pipeline['select2d_5'] = project.create_job(ws.uid,
                                        "select_2D",
                                        title = f"Select more than {thresh_particles} ptcl, better than {thresh_resolution} A resolution",
                                        params={
                                            "particle_count_above": thresh_particles,
                                            "resolution_better_than": thresh_resolution
                                            },
                                        connections={
                                            "particles": (pipeline['classify_good_template_picks2'].uid, "particles"),
                                            "templates": (pipeline['classify_good_template_picks2'].uid, "class_averages"),
                                            })
        pipeline['select2d_5'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['select2d_5'].uid} select2d_5")
        jobcounter += 1
        
    if resume == False or step == "classify_singleClass":
        resume = False
        
        pipeline['classify_singleClass'] = project.create_job(ws.uid,
                                                    "class_2D_new",
                                                    title="Group to single class",
                                                    connections={"particles": (pipeline['select2d_5'].uid, "particles_selected")},
                                                    params={
                                                        "class2D_K": 1,
                                                        "class2D_estimate_in_plane_pose": True,
                                                        "class2D_sigma_init_iter": 200,
                                                        "compute_use_ssd": scratch,
                                                        },
        )

        pipeline['classify_singleClass'].queue(lane=gpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['classify_singleClass'].uid} classify_singleClass")
        jobcounter += 1
        
    if resume == False or step == "filament_extrapolation":
        resume = False
        
        pipeline['classify_singleClass'].wait_for_done()
        particles = pipeline['classify_singleClass'].load_output('particles')

        resultgroup = pipeline['classify_singleClass'].download_file(f"{pipeline['classify_singleClass'].uid}_particles.csg")


        num_cpus = 1
        output_name = pipeline['classify_singleClass'].uid + "_extrapolated_filaments"
        input_particles = f"{pipeline['classify_singleClass'].uid}_020_particles.cs"
        passthrough = f"{pipeline['classify_singleClass'].uid}_passthrough_particles.cs"

        mics = particles.split_by("location/micrograph_path")

        # Initial JIT compilation before multiprocessing
        start = time.time()
        small_mic = max(mics.values(), key=len).slice(10)
        extrapolate_filaments(small_mic)
        tqdm.write(f"Finished numba pre-compiling ({time.time() - start:.4f}s)")

        tqdm.write("Processing remaining micrographs")

        with mp.Pool(min(len(mics), num_cpus)) as pool:
            extrapolated = list(
                tqdm(pool.imap_unordered(extrapolate_filaments, mics.values()), total=len(mics))
            )

        tqdm.write("Writing final results group")

        result = Dataset.append_many(*extrapolated)

        if not output_name:
            output_name = f"{input_particles[:-3]}_extrapolated"

        result.save(f"{output_name}.cs")
        with open(resultgroup) as f:
            csg = yaml.safe_load(f.read())
        new_results_group = {
            "location": csg["results"]["location"],
            "filament": csg["results"]["filament"],
        }
        new_results_group["location"]["metafile"] = new_results_group["location"][
            "metafile"
        ].replace(passthrough, f"{output_name}.cs")
        new_results_group["filament"]["metafile"] = new_results_group["filament"][
            "metafile"
        ].replace(input_particles, f"{output_name}.cs")
        new_csg = {
            "created": csg["created"],
            "group": csg["group"],
            "results": new_results_group,
            "version": csg["version"],
        }
        with open(f"{output_name}.csg", "w") as f:
            yaml.dump(new_csg, f)

        tqdm.write(f"Written {output_name}.csg")

        project.upload(f"exports/filamentExtrapolation/{pipeline['classify_singleClass'].uid}/{output_name}.cs", f'{output_name}.cs', overwrite=True)
        project.upload(f"exports/filamentExtrapolation/{pipeline['classify_singleClass'].uid}/{output_name}.csg", f'{output_name}.csg', overwrite=True)


        pipeline['importFilaments'] = project.create_job(ws.uid, "import_result_group",
                                        title="Import extrapolated filaments",
                                        params={
                                            "blob_path": os.path.abspath(f"{project.dir()}/exports/filamentExtrapolation/{pipeline['classify_singleClass'].uid}/{output_name}.csg"),
                                        }
                                        )
        pipeline['importFilaments'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['importFilaments'].uid} importFilaments")
        jobcounter += 1
        
    stage_apix = 2

    if resume == False or step == "extractFilaments":    
        resume = False
        
        pipeline['importFilaments'].wait_for_done()
        print(f"Stage 2, processing at {stage_apix} A/px")

        pipeline['extract4'] = project.create_job(ws.uid, 
                                    "extract_micrographs_cpu_parallel", 
                                    title="4th Extract and bin to 2A/px",
                                    connections={
                                        "micrographs": (pipeline['inspect1'].uid, "micrographs"),
                                        "particles": (pipeline['importFilaments'].uid, "particles")
                                    }, 
                                    params = {
                                        "box_size_pix" : find_box_size(apix_original, tube_diameter=tube_diameter, factor=particle_padding_factor),
                                        "bin_size_pix" : find_box_size(stage_apix, tube_diameter=tube_diameter, factor=particle_padding_factor),
                                    })
        pipeline['extract4'].queue(lane = cpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['extract4'].uid} extract4")
        jobcounter += 1
        
    if resume == False or step == "classify_extrapolated1":
        resume = False
        
        pipeline['extract4'].wait_for_done()
        pipeline['classify_extrapolated1'] = project.create_job(ws.uid,
                                                "class_2D_new",
                                                title="Extrapolated filaments classify 1",
                                                connections={"particles": (pipeline['extract4'].uid, "particles")},
                                                params={
                                                    "class2D_K": 50,
                                                    "class2D_estimate_in_plane_pose": True,
                                                    "class2D_num_full_iter": 2,
                                                    "class2D_sigma_init_iter": 200,
                                                    "compute_use_ssd": scratch,
                                                    },
        )
        pipeline['classify_extrapolated1'].queue(lane=gpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['classify_extrapolated1'].uid} classify_extrapolated1")
        jobcounter += 1
        
    if resume == False or step == "select2d_6":
        resume = False
        
        pipeline['classify_extrapolated1'].wait_for_done()
        
        particle_number = get_particle_number(pipeline['classify_extrapolated1'])
        if particle_number < 10000:
            nq_frac = 8
        else:
            nq_frac = 4
        
        thresh_particles, thresh_resolution = set_selection_limits(
            pipeline['classify_extrapolated1'], 
            dataset_fraction_per_class = 0.02, 
            nyquist_fraction = nq_frac)
        
        pipeline['select2d_6'] = project.create_job(ws.uid, 
                                        "select_2D",
                                        title = f"Select more than {thresh_particles} ptcl, better than {thresh_resolution} A resolution",
                                        params={
                                            "particle_count_above": thresh_particles,
                                            "resolution_better_than": thresh_resolution
                                            },
                                        connections={
                                            "particles": (pipeline['classify_extrapolated1'].uid, "particles"),
                                            "templates": (pipeline['classify_extrapolated1'].uid, "class_averages"),
                                            })
        pipeline['select2d_6'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['select2d_6'].uid} select2d_6")
        jobcounter += 1
        
    if resume == False or step == "classify_extrapolated2":
        resume = False
        pipeline['classify_extrapolated2'] = project.create_job(ws.uid,
                                                    "class_2D_new",
                                                    title="Extrapolated filaments classify 2",
                                                    connections={"particles": (pipeline['select2d_6'].uid, "particles_selected")},
                                                    params={
                                                        "class2D_K": 50,
                                                        "class2D_estimate_in_plane_pose": True,
                                                        "class2D_num_full_iter": 2,
                                                        "class2D_sigma_init_iter": 200,
                                                        "compute_use_ssd": scratch,
                                                        },
        )

        pipeline['classify_extrapolated2'].queue(lane = gpu_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['classify_extrapolated2'].uid} classify_extrapolated2")
        jobcounter += 1
        
    if resume == False or step == "select2d_7":
        resume = False
        
        pipeline['classify_extrapolated2'].wait_for_done()
        
        particle_number = get_particle_number(pipeline['classify_extrapolated2'])
        if particle_number < 10000:
            nq_frac = 8
        else:
            nq_frac = 4
        
        thresh_particles, thresh_resolution = set_selection_limits(
            pipeline['classify_extrapolated2'], 
            dataset_fraction_per_class = 0.02, 
            nyquist_fraction = nq_frac)
            
        pipeline['select2d_7'] = project.create_job(ws.uid, 
                                        "select_2D",
                                        title = f"Select more than {thresh_particles} ptcl, better than {thresh_resolution} resolution",
                                        params={
                                            "particle_count_above": thresh_particles,
                                            "resolution_better_than": thresh_resolution
                                            },
                                        connections={
                                            "particles": (pipeline['classify_extrapolated2'].uid, "particles"),
                                            "templates": (pipeline['classify_extrapolated2'].uid, "class_averages"),
                                            })
        pipeline['select2d_7'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['select2d_7'].uid} select2d_7")
        jobcounter += 1
        
    if resume == False or step == "models":
        resume = False
        
        pipeline['models'] = project.create_job(ws.uid, "import_volumes")

        pipeline['models'].set_param("volume_blob_path", models)
        pipeline['models'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['models'].uid} models")
        jobcounter += 1
        
    if resume == False or step == "protofilament_classification":
        resume = False
        pipeline['models'].wait_for_done()
        pipeline['protofilament_classification'] = project.create_job(ws.uid,
                                                    "hetero_refine",
                                                    title="Classify protofilament classes on imported maps",
                                                    params={
                                                        "multirefine_force_hard_class" : True,
                                                        "compute_use_ssd": scratch,
                                                        "multirefine_batch_size_per_class": 10000,
                                                        },
        )

        pipeline['protofilament_classification'].connect("particles", pipeline['select2d_7'].uid, "particles_selected")
        for i in range(6):
            pipeline['protofilament_classification'].connect("volume", pipeline['models'].uid, f"imported_volume_{i+1}")

        pipeline['protofilament_classification'].queue(lane = gpu_highmem_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['protofilament_classification'].uid} protofilament_classification")
        jobcounter += 1
        
    if resume == False or step == "importSortedFilaments":
        resume = False
        pipeline['protofilament_classification'].wait_for_done()
        
        particles = pipeline['protofilament_classification'].load_output("particles_all_classes")
        #particles_classified_csg = os.path.abspath(f"{project.dir()}/exports/groups/{pipeline['protofilament_classification'].uid}_particles_all_classes/{pipeline['protofilament_classification'].uid}_particles_all_classes_exported.csg")
        
        input_particles = f"{pipeline['protofilament_classification'].uid}_particles_all_classes.cs"
        output_name = f"{pipeline['protofilament_classification'].uid}_particles_all_classes_assigned"
        
        particles_classified_csg = pipeline['protofilament_classification'].download_file(f"{pipeline['protofilament_classification'].uid}_particles_all_classes.csg", f"{pipeline['protofilament_classification'].uid}_particles_all_classes.csg")
        
        conf_threshold = 0.7
        mics = particles.split_by("location/micrograph_path")

        new_tubes = []

        for mic_name in tqdm(mics):
            mic = mics[mic_name]

            fil_uids = np.unique(mic["filament/filament_uid"])
            new_fil_uids = fil_uids.tolist()

            for fil_uid in fil_uids:
                tube = mic.query({"filament/filament_uid": fil_uid})
                classes = class_posterior(tube)

                smoothened = []
                for i in range(0, len(tube)):
                    window = classes[
                        np.clip(i - 3, 0, None) : np.clip(i + 4, None, len(tube))
                    ]
                    mode, count = stats.mode(window)
                    smoothened.append(mode)

                changes = np.diff(smoothened)
                split_pos = (
                    [0] + (np.argwhere(changes) + 1).flatten().tolist() + [len(tube)]
                )

                for i in range(len(split_pos) - 1):
                    split = tube.take(range(split_pos[i], split_pos[i + 1]))
                    modal_class = smoothened[split_pos[i]]

                    if (
                        np.count_nonzero(
                            classes[split_pos[i] : split_pos[i + 1]] == modal_class
                        )
                        / len(split)
                        < conf_threshold
                    ):
                        continue

                    split["alignments3D_multi/class_posterior"] = (
                        split["alignments3D_multi/class"] == modal_class
                    ).astype(int)

                    new_fil_uid = max(new_fil_uids) + 1
                    new_fil_uids.append(new_fil_uid)  # added to prevent bug where 1 filament per micrograph was generated
                    
                    split["filament/filament_uid"] = new_fil_uid
                    new_tubes.append(split)

        new_tubes_cs = Dataset.append_many(*new_tubes)

        if not output_name:
            output_name = f"{input_particles[:-3]}_assigned"

        new_tubes_cs.save(f"{output_name}.cs")
        csg = []
        with open(particles_classified_csg) as f:
            for line in f:
                if "metafile:" not in line:
                    csg.append(line)
                else:
                    csg.append(f"    metafile: \'>{output_name}.cs\'\n")
        csg = ''.join(csg)
        #csg = csg.replace(f"{pipeline['protofilament_classification'].uid}_passthrough_particles_all_classes.cs", f"{output_name}.cs")
        with open(f"{output_name}.csg", "w") as f:
            f.write(csg)

        tqdm.write(f"Written {output_name}.csg")
        
        init_assignments = np.sum(particles['alignments3D_multi/class_posterior'], axis=0)
        final_assignments = np.sum(new_tubes_cs['alignments3D_multi/class_posterior'], axis=0)
        for i in range(len(init_assignments)):
            print(f"Initial class {i} with {i+11} protofilaments had {init_assignments[i]} particles, final class {i} has {final_assignments[i]} particles.")
        
        project.upload(f"exports/groups/{pipeline['protofilament_classification'].uid}_particles_all_classes/{output_name}.cs", f'{output_name}.cs', overwrite=True)
        project.upload(f"exports/groups/{pipeline['protofilament_classification'].uid}_particles_all_classes/{output_name}.csg", f'{output_name}.csg', overwrite=True)      
        
        pipeline['importSortedFilaments']=project.create_job(ws.uid, "import_result_group",
                                                            title="Import sorted protofilaments",
                                                            params={
                                                                "blob_path": os.path.abspath(f"{project.dir()}/exports/groups/{pipeline['protofilament_classification'].uid}_particles_all_classes/{output_name}.csg"),
                                                            }
                                                        )
        pipeline['importSortedFilaments'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['importSortedFilaments'].uid} importSortedFilaments")
        jobcounter += 1
    
    if resume == False or step == "reconstructClassified":
        resume = False
        
        pipeline['importSortedFilaments'].wait_for_done()
        pipeline['reconstructClassified'] = project.create_job(ws.uid, "hetero_reconstruct_new", 
                                                title="Reconstruct the pf classes",
                                                connections={
                                                    "particles": (pipeline['importSortedFilaments'].uid, "particles_all_classes")
                                                },
                                                params={
                                                    "compute_use_ssd": scratch,
                                                }
                                                )
        pipeline['reconstructClassified'].queue(lane = gpu_highmem_short_lane)
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['reconstructClassified'].uid} reconstructClassified")
    if resume == False or step == "splitVolumes":
        resume = False
        
        pipeline['splitVolumes'] = project.create_job(ws.uid, 
                                    "split_volumes_group",
                                            connections = {
                                            "particles" : (pipeline['reconstructClassified'].uid, "particles_all_classes"),
                                            "volumes_all_classes" : (pipeline['reconstructClassified'].uid, "volumes_all_classes")
                                    })
        pipeline['splitVolumes'].queue()
        log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline['splitVolumes'].uid} splitVolumes")
        jobcounter += 1
        
    if resume == False or step.startswith( "helical_refine"):
        resume = False
        
        pipeline['splitVolumes'].wait_for_done()
        
        valid_classes = []

        for i, pf_start in enumerate(microtubule_params_theory):
            
            print(f"Processing {i+11} protofilament class, rise of {microtubule_params_theory[pf_start]['rise']} and twist of {microtubule_params_theory[pf_start]['twist']}")
            class_size = len(pipeline['splitVolumes'].load_output(f'particles_class_{i}'))
            print(f"Class has {class_size} particles.")
            if class_size < 5000:
                print(f"Class {i+11} has less than 5000 particles, skipping.")
                continue
            else:
                pipeline[f'reconstruct_{pf_start}_protofilaments'] = project.create_job(ws.uid,
                                                                                    "homo_reconstruct",
                                                                                    title = f"Homogeneous reconstruction for {pf_start} protofilament class",
                                                                                    connections = {
                                                                                        "particles" : (pipeline['splitVolumes'].uid, f"particles_class_{i}"),
                                                                                    },
                                                                                    params = {
                                                                                        "compute_use_ssd": scratch,
                                                                                    }
                                                                                )
                pipeline[f'reconstruct_{pf_start}_protofilaments'].queue(lane   = gpu_highmem_short_lane)
                                                                        
                pipeline[f'helical_refine_{pf_start}'] = project.create_job(ws.uid, "helix_refine",
                                                    title=f"Helical refinement for {pf_start} protofilament class",
                                                    params={
                                                        "nu_refine": True,
                                                        "refine_scale_min": True,
                                                        "refine_hsym_order" : i+11,
                                                        "refine_init_shift": microtubule_params_theory[pf_start]['rise'],
                                                        "refine_init_twist": microtubule_params_theory[pf_start]['twist'],
                                                        "compute_use_ssd": scratch,
                                                    },
                                                    connections = {
                                                        "particles" : (pipeline['splitVolumes'].uid, f"particles_class_{i}"),
                                                        "volume" : (pipeline[f'reconstruct_{pf_start}_protofilaments'].uid, f"volume")
                                                    })
                pipeline[f'helical_refine_{pf_start}'].queue(lane = gpu_highmem_lane)
                log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'helical_refine_{pf_start}'].uid} helical_refine_{pf_start}")
                jobcounter += 1
                valid_classes.append(pf_start)
                
        log_append(logfile_path, f"Valid classes are: {[i for i in valid_classes]}")
    
    stage_apix = 1.3
    print(f"Processing at {stage_apix} A/px")
    
    if resume == False or step.startswith("extract5"):
        resume = False
        for pf_start in valid_classes:
            try:
                #pipeline[f'helical_refine_{i+11}'].wait_for_done()
                pipeline[f'extract5_{pf_start}_pf'] = project.create_job(ws.uid,
                                                "extract_micrographs_cpu_parallel",
                                                title=f"Extract at {stage_apix} A/px for {pf_start} protofilament class",
                                                connections={
                                                    "micrographs": (pipeline[f'inspect1'].uid, "micrographs"),
                                                    "particles": (pipeline[f'helical_refine_{pf_start}'].uid, "particles")
                                                },
                                                params={
                                                    "box_size_pix": find_box_size(apix_original, particle_padding_factor, tube_diameter),
                                                    "bin_size_pix": find_box_size(stage_apix, particle_padding_factor, tube_diameter),
                                                    "recenter_using_shifts": True,
                                                    "output_f16": True,
                                                    "compute_num_cores": 24
                                                })
                pipeline[f'extract5_{pf_start}_pf'].queue(lane = cpu_lane)
                log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'extract5_{pf_start}_pf'].uid} extract5_{pf_start}_pf")
                jobcounter += 1 
                
            except KeyError:
                print(f"Class {pf_start} discarded")
                continue
    
    if resume == False or step.startswith("fullsize_helical_refine"):
        resume = False
        for pf_start in valid_classes:
            pf = int(pf_start.split("-")[0])
            #pipeline[f'extract5_{pf}_pf'].wait_for_done()
            pipeline[f'fullsize_helical_refine_{pf_start}'] = project.create_job(ws.uid, "helix_refine",
                                                    title=f"Helical refinement for {pf_start} protofilament class at full size",
                                                    params={
                                                        "nu_refine": True,
                                                        "refine_scale_min": True,
                                                        "refine_hsym_order" : pf,
                                                        "refine_init_shift": microtubule_params_theory[pf_start]['rise'],
                                                        "refine_init_twist": microtubule_params_theory[pf_start]['twist'],
                                                        "compute_use_ssd": scratch,
                                                        "low_memory_mode": True
                                                    },
                                                    connections = {
                                                        "particles" : (pipeline[f'extract5_{pf_start}_pf'].uid, "particles"),
                                                        "volume" : (pipeline[f'helical_refine_{pf_start}'].uid, "volume")
                                                    })
            pipeline[f'fullsize_helical_refine_{pf_start}'].queue(lane = gpu_veryhighmem_lane)
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'fullsize_helical_refine_{pf_start}'].uid} fullsize_helical_refine_{pf_start}")
            jobcounter += 1
            
    if resume == False or step.startswith("ctf_refine_per_particle"):
        resume = False
        for pf_start in valid_classes:
            #pf, start = i.split("-")
            #pipeline[f'fullsize_helical_refine_{pf_start}'].wait_for_done()
            pipeline[f'ctf_refine_per_particle_{pf_start}'] = project.create_job(ws.uid, "ctf_refine_local",
                                                    title=f"CTF refinement for {pf_start} protofilament class",
                                                    connections = {
                                                        "particles" : (pipeline[f'fullsize_helical_refine_{pf_start}'].uid, "particles"),
                                                        "volume" : (pipeline[f'fullsize_helical_refine_{pf_start}'].uid, "volume")
                                                    },
                                                    params = {
                                                        "compute_use_ssd": scratch,
                                                    })
            pipeline[f'ctf_refine_per_particle_{pf_start}'].queue(lane = gpu_highmem_short_lane)
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'ctf_refine_per_particle_{pf_start}'].uid} ctf_refine_per_particle_{pf_start}")
            jobcounter += 1
    
    if resume == False or step.startswith("fullsize_c1_helix_refine"):
        resume = False
        
        for pf_start in valid_classes:
            pf = int(pf_start.split("-")[0])
            #pipeline[f'extract5_{pf}_pf'].wait_for_done()
            pipeline[f'fullsize_c1_helix_refine_{pf_start}'] = project.create_job(ws.uid, "helix_refine",
                                                    title=f"C1 refinement for {pf_start} protofilament class at full size",
                                                    params={
                                                        "nu_refine": True,
                                                        "refine_scale_min": True,
                                                        "refine_hsym_order" : 1,
                                                        "refine_init_shift": microtubule_params_theory[pf_start]['rise'],
                                                        "refine_init_twist": microtubule_params_theory[pf_start]['twist'],
                                                        "compute_use_ssd": scratch,
                                                        "low_memory_mode": True 
                                                    },
                                                    connections = {
                                                        "particles" : (pipeline[f'ctf_refine_per_particle_{pf_start}'].uid, "particles"),
                                                        "volume" : (pipeline[f'fullsize_helical_refine_{pf_start}'].uid, "volume")
                                                    })
            pipeline[f'fullsize_c1_helix_refine_{pf_start}'].queue(lane = gpu_veryhighmem_lane)
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid} fullsize_c1_helix_refine_{pf_start}")
            jobcounter += 1
            
    if resume == False or step.startswith("psi_unify"):
        resume = False
        
        for pf_start in valid_classes:
            pipeline[f'fullsize_c1_helix_refine_{pf_start}'].wait_for_done()
            particles = pipeline[f'fullsize_c1_helix_refine_{pf_start}'].load_output("particles")
            input_particles = f"{pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid}_particles.cs"
            output_name = f"{pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid}_particles_psi_corrected"
            particles_classified_csg = pipeline[f'fullsize_c1_helix_refine_{pf_start}'].download_file(f"{pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid}_particles.csg", f"{pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid}_particles.csg")

                        
            mics = particles.split_by("location/micrograph_path")

            unified_tubes = []

            for mic_name in tqdm(mics):
                tubes = mics[mic_name].split_by("filament/filament_uid")

                for filament_uid in tubes:
                    tube = tubes[filament_uid]

                    if len(tube) < 3:
                        continue

                    eulers = R.from_rotvec(tube["alignments3D/pose"]).as_euler(
                        "ZYZ", degrees=True
                    )
                    psi = eulers[:, 2]

                    if not (abs(np.diff(psi)) > 1).any():
                        unified_tubes.append(tube)
                        continue

                    vals, bins = np.histogram(psi, bins="doane", range=(-180, 180))
                    mode_ind = vals.argmax()

                    modal_mask = (psi > bins[mode_ind] - 10) & (psi < bins[mode_ind + 1] + 10)
                    modal_points = psi[modal_mask]

                    if len(modal_points) < 3:
                        continue

                    fit = (
                        np.polynomial.polynomial.Polynomial.fit(
                            np.argwhere(modal_mask).flatten(),
                            modal_points,
                            2,
                        )
                        .convert()
                        .coef
                    )
                    corrected_psi = np.polynomial.polynomial.polyval(np.arange(len(tube)), fit)
                    eulers[:, 2] = corrected_psi
                    poses = R.from_euler("ZYZ", eulers, degrees=True).as_rotvec()
                    tube["alignments3D/pose"] = poses

                    unified_tubes.append(tube)

            unified_tubes_cs = Dataset.append_many(*unified_tubes)

            if not output_name:
                output_name = f"{input_particles[:-3]}_upsi"

            unified_tubes_cs.save(f"{output_name}.cs")

            csg = []
            with open(input_particles + "g") as f:
                for line in f:
                    if "metafile:" not in line:
                        csg.append(line)
                    else:
                        csg.append(f"    metafile: \'>{output_name}.cs\'\n")
            csg = ''.join(csg)
            #csg = csg.replace(input_particles, f"{output_name}.cs")
            with open(f"{output_name}.csg", "w") as f:
                f.write(csg)


            tqdm.write(f"Written {output_name}.csg")

            project.upload(f"exports/groups/{pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid}_particles/{output_name}.cs", f'{output_name}.cs', overwrite=True)
            project.upload(f"exports/groups/{pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid}_particles/{output_name}.csg", f'{output_name}.csg', overwrite=True)  
            
            pipeline[f'psi_unify_{pf_start}']=project.create_job(ws.uid, "import_result_group",
                                   title=f"Import psi-unified particles for {pf_start} protofilament class",
                                   params={
                                       "blob_path": os.path.abspath(f"{project.dir()}/exports/groups/{pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid}_particles/{output_name}.csg"),
                                   }
                                   )
            pipeline[f'psi_unify_{pf_start}'].queue()
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'psi_unify_{pf_start}'].uid} psi_unify_{pf_start}")
            jobcounter += 1
            
    if resume == False or step.startswith("local_refine_psi"):
        resume = False
        
        for pf_start in valid_classes:
            pipeline[f'psi_unify_{pf_start}'].wait_for_done()
            pipeline[f'local_refine_psi_{pf_start}'] = project.create_job(ws.uid, "new_local_refine",
                                                    title=f"Local refinement for {pf_start} psi-unified protofilament class",
                                                    params={
                                                        "use_alignment_prior": True,
                                                        "reinitialize_rs": True,
                                                        "reinitialize_ss": True,
                                                        "sigma_prior_r": 5,
                                                        "sigma_prior_s": 3,
                                                        "compute_use_ssd": scratch,
                                                        #"low_memory_mode": True 
                                                    },
                                                    connections = {
                                                        "particles" : (pipeline[f'psi_unify_{pf_start}'].uid, "particles"),
                                                        "volume" : (pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid, "volume"),
                                                        "mask": (pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid, "mask")
                                                    })
            pipeline[f'local_refine_psi_{pf_start}'].queue(lane = gpu_highmem_lane)
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'local_refine_psi_{pf_start}'].uid} local_refine_psi_{pf_start}")
            jobcounter += 1
    
    if resume == False or step.startswith("unify_phi_1"):
        resume = False
        
        for pf_start in valid_classes:
            pipeline[f'local_refine_psi_{pf_start}'].wait_for_done()
            particles = pipeline[f'local_refine_psi_{pf_start}'].load_output("particles")
            input_particles = f"{pipeline[f'local_refine_psi_{pf_start}'].uid}_particles.cs"
            output_name = f"{pipeline[f'local_refine_psi_{pf_start}'].uid}_particles_phi_corrected1"
            particles_classified_csg = pipeline[f'local_refine_psi_{pf_start}'].download_file(f"{pipeline[f'local_refine_psi_{pf_start}'].uid}_particles.csg", f"{pipeline[f'local_refine_psi_{pf_start}'].uid}_particles.csg")

            mics = particles.split_by("location/micrograph_path")

            unified_tubes = []

            pfn, _ = pf_start.split("-")
            pfn = int(pfn)
            
            num_cpus = 6

            # Here added from new version N+17 lines
            pdf = PdfPages(f"{output_name}.pdf")

            with mp.Pool(min(len(mics), num_cpus)) as pool:
                for figs, tube in tqdm(
                    pool.imap_unordered(
                        phi_process_mic,
                        [(mic, pfn) for mic in mics.values()],
                        chunksize=1,
                    ),
                    total=len(mics),
                ):
                    unified_tubes.append(tube)
                    for fig in figs:
                        pdf.savefig(fig)
                        plt.close(fig)

            pdf.close()
            unified_tubes_cs = Dataset.append_many(*unified_tubes)

            if not output_name:
                output_name = f"{input_particles[:-3]}_uphi"

            unified_tubes_cs.save(f"{output_name}.cs")

            csg = []
            with open(input_particles + "g") as f:
                for line in f:
                    if "metafile:" not in line:
                        csg.append(line)
                    else:
                        csg.append(f"    metafile: \'>{output_name}.cs\'\n")
            csg = ''.join(csg)
            #csg = csg.replace(input_particles, f"{output_name}.cs")
            with open(f"{output_name}.csg", "w") as f:
                f.write(csg)

            tqdm.write(f"Written {output_name}.csg")
            
            project.upload(f"exports/groups/{pipeline[f'local_refine_psi_{pf_start}'].uid}_particles/{output_name}.cs", f'{output_name}.cs', overwrite=True)
            project.upload(f"exports/groups/{pipeline[f'local_refine_psi_{pf_start}'].uid}_particles/{output_name}.csg", f'{output_name}.csg', overwrite=True)
            
            pipeline[f'unify_phi_1_{pf_start}']=project.create_job(ws.uid, "import_result_group",
                                   title=f"Import phi-unified particles for {pf_start} protofilament class",
                                   params={
                                       "blob_path": os.path.abspath(f"{project.dir()}/exports/groups/{pipeline[f'local_refine_psi_{pf_start}'].uid}_particles/{output_name}.csg"),
                                   }
                                   )
            pipeline[f'unify_phi_1_{pf_start}'].queue()
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'unify_phi_1_{pf_start}'].uid} unify_phi_1_{pf_start}")
            jobcounter += 1
            
    if resume == False or step.startswith("local_refine_phi_1"):
        resume = False
        
        for pf_start in valid_classes:
            pipeline[f'unify_phi_1_{pf_start}'].wait_for_done()
            pipeline[f'local_refine_phi_1_{pf_start}'] = project.create_job(ws.uid, "new_local_refine",
                                                    title=f"Local refinement for {pf_start} phi-unified protofilament class",
                                                    params={
                                                        "use_alignment_prior": True,
                                                        "reinitialize_rs": True,
                                                        "reinitialize_ss": True,
                                                        "sigma_prior_r": 5,
                                                        "sigma_prior_s": 3,
                                                        "compute_use_ssd": scratch,
                                                        #"low_memory_mode": True 
                                                    },
                                                    connections = {
                                                        "particles" : (pipeline[f'unify_phi_1_{pf_start}'].uid, "particles"),
                                                        "volume" : (pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid, "volume"),
                                                        "mask": (pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid, "mask")
                                                    })
            pipeline[f'local_refine_phi_1_{pf_start}'].queue(lane = gpu_highmem_lane)
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'local_refine_phi_1_{pf_start}'].uid} local_refine_phi_1_{pf_start}")
            jobcounter += 1
            
    if resume == False or step.startswith("unify_phi_2"):
        resume = False
        
        for pf_start in valid_classes:
            pipeline[f'local_refine_phi_1_{pf_start}'].wait_for_done()
            particles = pipeline[f'local_refine_phi_1_{pf_start}'].load_output("particles")
            input_particles = f"{pipeline[f'local_refine_phi_1_{pf_start}'].uid}_particles.cs"
            output_name = f"{pipeline[f'local_refine_phi_1_{pf_start}'].uid}_particles_phi_corrected2"
            particles_classified_csg = pipeline[f'local_refine_phi_1_{pf_start}'].download_file(f"{pipeline[f'local_refine_phi_1_{pf_start}'].uid}_particles.csg", f"{pipeline[f'local_refine_phi_1_{pf_start}'].uid}_particles.csg")

            mics = particles.split_by("location/micrograph_path")

            unified_tubes = []

            pfn, _ = pf_start.split("-")
            pfn = int(pfn)
            num_cpus = 6

            pdf = PdfPages(f"{output_name}.pdf")

            with mp.Pool(min(len(mics), num_cpus)) as pool:
                for figs, tube in tqdm(
                    pool.imap_unordered(
                        phi_process_mic,
                        [(mic, pfn) for mic in mics.values()],
                        chunksize=1,
                    ),
                    total=len(mics),
                ):
                    unified_tubes.append(tube)
                    for fig in figs:
                        pdf.savefig(fig)
                        plt.close(fig)

            pdf.close()
            unified_tubes_cs = Dataset.append_many(*unified_tubes)

            if not output_name:
                output_name = f"{input_particles[:-3]}_uphi"

            unified_tubes_cs.save(f"{output_name}.cs")
            
            csg = []
            with open(input_particles + "g") as f:
                for line in f:
                    if "metafile:" not in line:
                        csg.append(line)
                    else:
                        csg.append(f"    metafile: \'>{output_name}.cs\'\n")
            csg = ''.join(csg)
            #csg = csg.replace(input_particles, f"{output_name}.cs")
            with open(f"{output_name}.csg", "w") as f:
                f.write(csg)

            tqdm.write(f"Written {output_name}.csg")
            
            project.upload(f"exports/groups/{pipeline[f'local_refine_phi_1_{pf_start}'].uid}_particles/{output_name}.cs", f'{output_name}.cs', overwrite=True)
            project.upload(f"exports/groups/{pipeline[f'local_refine_phi_1_{pf_start}'].uid}_particles/{output_name}.csg", f'{output_name}.csg', overwrite=True)
            
            pipeline[f'unify_phi_2_{pf_start}']=project.create_job(ws.uid, "import_result_group",
                                   title=f"Import phi-unified particles for {pf_start} protofilament class",
                                   params={
                                       "blob_path": os.path.abspath(f"{project.dir()}/exports/groups/{pipeline[f'local_refine_phi_1_{pf_start}'].uid}_particles/{output_name}.csg"),
                                   }
                                   )
            pipeline[f'unify_phi_2_{pf_start}'].queue()
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'unify_phi_2_{pf_start}'].uid} unify_phi_2_{pf_start}")
            jobcounter += 1
            
    if resume == False or step.startswith("local_refine_phi_2"):
        resume = False
        
        for pf_start in valid_classes:
            pipeline[f'unify_phi_2_{pf_start}'].wait_for_done()
            pipeline[f'local_refine_phi_2_{pf_start}'] = project.create_job(ws.uid, "new_local_refine",
                                                    title=f"Local refinement for {pf_start} phi-unified protofilament class",
                                                    params={
                                                        "use_alignment_prior": True,
                                                        "sigma_prior_r": 3,
                                                        "sigma_prior_s": 2,
                                                        "compute_use_ssd": scratch,
                                                        #"low_memory_mode": True 
                                                    },
                                                    connections = {
                                                        "particles" : (pipeline[f'unify_phi_2_{pf_start}'].uid, "particles"),
                                                        "volume" : (pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid, "volume"),
                                                        "mask": (pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid, "mask")
                                                    })
            pipeline[f'local_refine_phi_2_{pf_start}'].queue(lane = gpu_highmem_lane)
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'local_refine_phi_2_{pf_start}'].uid} local_refine_phi_2_{pf_start}")
            jobcounter += 1
            
    if resume == False or step.startswith("symmetry_expansion"):
        resume = False
        for pf_start in valid_classes:
            pipeline[f'local_refine_phi_2_{pf_start}'].wait_for_done()
            pf, start = pf_start.split("-")
            order, rise, twist = get_helical_refinement_results(pipeline[f'fullsize_helical_refine_{pf_start}'])
            pipeline[f'symmetry_expansion_{pf_start}'] = project.create_job(ws.uid, "sym_expand",
                                                    title=f"Symmetry expansion for {pf_start} phi-unified protofilament class",
                                                    params={
                                                        "sym_num_rises": order,
                                                        "sym_rise_A": rise,
                                                        "sym_twist_deg": twist
                                                    },
                                                    connections = {
                                                        "particles" : (pipeline[f'local_refine_phi_2_{pf_start}'].uid, "particles"),
                                                    })
            pipeline[f'symmetry_expansion_{pf_start}'].queue(lane=cpu_highmem_lane)
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'symmetry_expansion_{pf_start}'].uid} symmetry_expansion_{pf_start}")
            
    if resume == False or step.startswith("local_refine_symmetry_expansion"):
        resume = False
        for pf_start in valid_classes:
            pipeline[f'symmetry_expansion_{pf_start}'].wait_for_done()
            pipeline[f'local_refine_symmetry_expansion_{pf_start}'] = project.create_job(ws.uid, "new_local_refine",
                                                    title=f"Local refinement for {pf_start} symmetry-expanded protofilament class",
                                                    params={
                                                        "use_alignment_prior": True,
                                                        "sigma_prior_r": 3,
                                                        "sigma_prior_s": 2,
                                                        "compute_use_ssd": scratch,
                                                        #"low_memory_mode": True 
                                                    },
                                                    connections = {
                                                        "particles" : (pipeline[f'symmetry_expansion_{pf_start}'].uid, "particles"),
                                                        "volume" : (pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid, "volume"),
                                                        "mask": (pipeline[f'fullsize_c1_helix_refine_{pf_start}'].uid, "mask")
                                                    })
            pipeline[f'local_refine_symmetry_expansion_{pf_start}'].queue(lane=gpu_highmem_lane)
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'local_refine_symmetry_expansion_{pf_start}'].uid} local_refine_symmetry_expansion_{pf_start}")
            jobcounter += 1
            
    if resume == False or step.startswith("summary"):
        resume = False
        for pf_start in valid_classes:
            pipeline[f'local_refine_symmetry_expansion_{pf_start}'].wait_for_done()
            print(f"Processed {pf_start} protofilament class")
            print(f'Final resolution achieved: {get_refinement_resolution(pipeline[f"local_refine_symmetry_expansion_{pf_start}"])}')
            log_append(logfile_path, f"Queued {jobcounter:02d} {project.uid} {ws.uid} {pipeline[f'local_refine_symmetry_expansion_{pf_start}'].uid} summary_{pf_start}")
            jobcounter += 1
            