from cryosparc.tools import CryoSPARC
import yaml
import os

config_path = os.path.join(os.path.dirname(__file__), "cs_config.yml")
with open(config_path, "r") as f:
    cs_conf = yaml.safe_load(f)

cs = CryoSPARC(
    license=cs_conf["license"],
    host=cs_conf["host"],
    base_port=cs_conf["base_port"],
    email=cs_conf["email"],
    password=cs_conf["password"]
)

print("Testing CryoSPARC connection...")
assert cs.test_connection(), "Connection failed"
print("✅ Connection successful!")

