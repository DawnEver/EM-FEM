import os
import logging
import datetime

__all__ = ["log"]

current_time = datetime.datetime.now()
today = current_time.strftime("%Y-%m-%d")
log_folder_path = "./log/"
if not os.path.exists(log_folder_path):
    os.mkdir(log_folder_path)
today_log_filename = f"{log_folder_path}{today}.log"
log_file_format = "%(asctime)s %(levelname)s %(message)s \n\tlocation: %(filename)s %(funcName)s line%(lineno)d \n\tprocess: %(processName)s pid:%(process)d \n\tthread: %(threadName)s id:%(thread)d"
log_console_format = "%(message)s %(asctime)s"
log_date_format = "%m/%d/%Y %H:%M:%S"

logger = logging.getLogger()
logger.setLevel(logging.INFO)

logging_file_handler = logging.FileHandler(today_log_filename, encoding="utf-8")
logging_file_handler.setLevel(logging.INFO)
logging_file_format = logging.Formatter(log_file_format, datefmt=log_date_format)
logging_file_handler.setFormatter(logging_file_format)
logger.addHandler(logging_file_handler)


logging_console_handler = logging.StreamHandler()
logging_console_handler.setLevel(logging.INFO)
logging_console_format = logging.Formatter(log_console_format, datefmt=log_date_format)
logging_console_handler.setFormatter(logging_console_format)
logger.addHandler(logging_console_handler)


def log(msg: str, level: str = "INFO"):
    level = level.upper()
    match level:
        case "DEBUG":
            logging.debug(msg)
        case "INFO":
            logging.info(msg)
        case "WARNING":
            logging.warning(msg)
        case "ERROR":
            logging.error(msg)
        case "CRITICAL":
            logging.critical(msg)
        case _:
            logging.error("Unknown log level!")
