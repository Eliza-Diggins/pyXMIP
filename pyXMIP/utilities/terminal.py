import multiprocessing
import sys
import time

from pyXMIP.utilities.core import xsparams


class Spinner:
    busy = False
    delay = 0.05
    process = None
    stream = (
        sys.stdout
        if xsparams["system"]["logging"]["main"]["stream"] in ["STDOUT", "stdout"]
        else sys.stderr
    )
    formatter = xsparams["system"]["logging"]["main"]["format"]

    @staticmethod
    def spinning_cursor():
        while 1:
            for cursor in ["*", "**", "***", "****"]:
                yield cursor

    def __init__(self, text="", delay=None):
        self.text = text
        self.spinner_generator = self.spinning_cursor()
        if delay and float(delay):
            self.delay = delay

    def spinner_task(self):
        while self.busy:
            self.stream.write(
                self.formatter
                % dict(
                    message=f"{self.text} : {next(self.spinner_generator)}",
                    name="pyXs",
                    levelname="OPERATION",
                    asctime=time.asctime(),
                )
            )
            self.stream.flush()
            time.sleep(self.delay)
            self.stream.write("\r")
            self.stream.flush()

    def __enter__(self):
        self.busy = True
        self.process = multiprocessing.Process(target=self.spinner_task, daemon=True)
        self.process.start()

    def __exit__(self, exception, value, tb):
        self.process.kill()
        self.stream.write("\r")
        self.stream.write(
            self.formatter
            % dict(
                message=f"{self.text} : [COMPLETE]",
                name="pyXs",
                levelname="OPERATION",
                asctime=time.asctime(),
            )
        )
        self.stream.flush()
        self.busy = False
        time.sleep(self.delay)
        if exception is not None:
            return False
