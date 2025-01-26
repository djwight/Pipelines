import os
import logging


def init_logger(outdir: str) -> None:
    # logging formats
    log_fmt = "%(asctime)s [%(process)d] [%(levelname)s] %(message)s"
    date_fmt = "%Y-%m-%d %H:%M:%S"
    formatter = logging.Formatter(fmt=log_fmt, datefmt=date_fmt)

    # logging config
    logging.basicConfig(
        format=log_fmt,
        datefmt=date_fmt,
        level=logging.DEBUG,
        filename=os.path.join(outdir, "run.log"),
    )
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logging.getLogger().addHandler(stream_handler)
