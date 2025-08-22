import yaml
import subprocess

with open("config.yaml") as f:
    cfg = yaml.safe_load(f)

cmd = ["python", "exe_ecc.py"] + cfg["args"]
subprocess.run(cmd)
