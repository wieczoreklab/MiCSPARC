from cryosparc.tools import CryoSPARC
import os
import yaml
import argparse


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


parser = argparse.ArgumentParser()
parser.add_argument("--project", help="Project name")
parser.add_argument("--job", help="Job name")
parser.add_argument("--target_dir", required=False, help="Target directory")
args = parser.parse_args()

project = args.project
job = args.job
target_dir = args.target_dir


c3d = cs.find_job(project, job)

if target_dir:
    target_dir = args.target_dir
else:
    target_dir = f"{homedir}/Downloads/{project}_{c3d.uid}"

os.makedirs(f"{target_dir}", exist_ok=True)
for each in c3d.list_files():
    if each.endswith('_volume.mrc'):
        c3d.download_file(each, f"{target_dir}/{each}")