import logging
logger = logging.getLogger("main")
import time


_last_log_time = time.perf_counter()

def log(text, notime=False):
    global _last_log_time
    t = time.perf_counter()
    dt = t - _last_log_time
    _last_log_time = t
    if notime:
        logger.debug(f"{text}")
    else:
        logger.debug(f"{dt:4.2f} s | {text:60.60s}")
   
    