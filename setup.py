from setuptools import setup, find_packages, Distribution
from setuptools.command.build_py import build_py
import subprocess
import shutil
import os
from pathlib import Path


PACKAGE_NAME = 'scistree2'
# PROJECT_ROOT = str(Path(__file__).parent.resolve())
CPP_DIR = './src'
PACKAGE_DIR = 'scistree2'
# CPP_DIR = f'{PROJECT_ROOT}/src'
# PACKAGE_DIR = f'{PROJECT_ROOT}/scistree2'


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True  # <--- The magic line. Forces "platlib" (platform-specific) wheel.
    

class BuildWithMake(build_py):
    def run(self):
        self.build_cpp_with_make()
        super().run()
    
    def build_cpp_with_make(self):
        print("Building C++ project with make...")
        try:
            subprocess.run(['make'], cwd=CPP_DIR, check=True)
        except subprocess.CalledProcessError as e:
            print("Make failed.")
            raise e
        
        bin_dir = os.path.join(PACKAGE_DIR, 'bin')
        os.makedirs(bin_dir, exist_ok=True)
        shutil.copy(
            os.path.join(CPP_DIR, 'scistree'),
            os.path.join(bin_dir, 'scistree')
        )


setup(
    name=PACKAGE_NAME,
    packages=[PACKAGE_DIR],
    # cmdclass={'build_py': BuildWithMake},
    install_requires=[
        'numpy',
        'pptree',
        'phytreeviz'
    ],
    package_data={
        PACKAGE_NAME: ['bin/*'],  # Corrected path
    },
    include_package_data=False,
    distclass=BinaryDistribution,
)
