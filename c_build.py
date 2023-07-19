from setuptools import Extension
from setuptools.command.build_py import build_py as _build_py


class build_py(_build_py):
    def run(self):
        self.run_command("build_ext")
        return super().run()

    def initialize_options(self):
        super().initialize_options()
        if self.distribution.ext_modules == None:
            self.distribution.ext_modules = []

        self.distribution.ext_modules.append(
            Extension(
                "chromTools.c_io",
                sources=["chromTools/c_io.pyx"],
            )
        )