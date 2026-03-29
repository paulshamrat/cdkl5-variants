import os
import subprocess
from src.data_utils import load_config

class PalmettoBridge:
    def __init__(self, config=None):
        self.config = config if config else load_config()
        self.hpc_cfg = self.config.get('hpc', {}).get('palmetto', {})
        self.host = self.hpc_cfg.get('hostname')
        self.user = self.hpc_cfg.get('username')
        self.project_dir = self.hpc_cfg.get('project_dir')
        self.control_path = self.hpc_cfg.get('control_path')

    def _get_ssh_cmd(self, remote_cmd):
        cmd = ["ssh"]
        if self.control_path:
            cmd.extend(["-o", f"ControlPath={self.control_path}", "-o", "ControlMaster=no"])
        cmd.extend([f"{self.user}@{self.host}", remote_cmd])
        return cmd

    def list_remote_dir(self, path=None):
        """Lists files in a directory on Palmetto."""
        target_path = path if path else self.project_dir
        cmd = self._get_ssh_cmd(f"ls -F {target_path}")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return result.stdout.splitlines()
        except subprocess.CalledProcessError as e:
            return f"Error connecting to Palmetto: {e.stderr}"

    def generate_pbs_script(self, job_name, commands, output_path):
        """Generates a PBS script for Palmetto submission."""
        alloc = self.hpc_cfg.get('allocation', {})
        
        script_lines = [
            "#!/bin/bash",
            f"#PBS -N {job_name}",
            f"#PBS -l select={alloc.get('nodes')}:ncpus={alloc.get('cpus')}:mem={alloc.get('mem')}",
            f"#PBS -l walltime={alloc.get('walltime')}",
            "#PBS -j oe",
            f"cd {self.project_dir}",
            "\n"
        ]
        script_lines.extend(commands)
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            f.write("\n".join(script_lines))
        return output_path

    def run_remote_command(self, command):
        """Executes a arbitrary command on Palmetto via SSH."""
        cmd = self._get_ssh_cmd(command)
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return result.stdout
        except subprocess.CalledProcessError as e:
            return f"Error: {e.stderr}"

if __name__ == "__main__":
    # Quick test if run directly
    bridge = PalmettoBridge()
    print(f"Palmetto Project Directory: {bridge.project_dir}")
    print("Run `bridge.list_remote_dir()` to see remote files (requires SSH keys).")
