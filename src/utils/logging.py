import sys
from datetime import datetime
from pathlib import Path


def log_message(log_file: Path, message: str, level: str = "INFO") -> None:
    """Log a message with timestamp and level"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted_message = f"[{timestamp}] {level}: {message}"
    with open(log_file, "a") as log:
        log.write(formatted_message + "\n")

def log_script_command(info_file: Path) -> None:
    """Log the exact command used to execute the script"""
    command = f"python3 {sys.argv[0]} {' '.join(sys.argv[1:])}"
    log_message(info_file, f"Script command: {command}")
