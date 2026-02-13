import os
import sys
import logging
from datetime import datetime

def setup_logging(log_dir, log_name, log_level=logging.INFO):
    """
    Setup logging configuration to write to both file and console.
    
    Args:
        log_dir: Directory to save log file
        log_name: Name of the log file
        log_level: Logging level
    """
    os.makedirs(log_dir, exist_ok=True)
    
    # Create log filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"{log_name}_{timestamp}.log")
    
    # Create a logger instance
    logger = logging.getLogger(log_name)
    logger.setLevel(log_level)
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Create file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    
    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return log_file, logger