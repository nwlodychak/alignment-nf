import subprocess
from datetime import datetime
import logging

def run_command(command:str):
    """
    Execute a shell command and return process.
    :param command - shell command to run
    """
    print(f"Running: {command}")
    process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE, text=True)
    
    if process.returncode != 0:
        print(f"Error executing command: {command}")
        print(f"STDERR: {process.stderr}")
        return False
    return True