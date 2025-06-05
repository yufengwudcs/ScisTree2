from setuptools import setup
from setuptools.command.build_py import build_py
import subprocess
import shutil
import os

CPP_DIR = './src'
PACKAGE_NAME = 'scistree2'


class BuildWithCMake(build_py):
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
        
        bin_dir = os.path.join(PACKAGE_NAME, 'bin')
        os.makedirs(bin_dir, exist_ok=True)
        shutil.copy(
            os.path.join(CPP_DIR, 'scistree'),
            os.path.join(bin_dir, 'scistree')
        )


setup(
    name=PACKAGE_NAME,
    version='0.1.0',
    packages=[PACKAGE_NAME],
    cmdclass={'build_py': BuildWithCMake},
    install_requires=[
        'numpy',
        'pptree'
    ],
    python_requires='>=3.6',
    package_data={
        PACKAGE_NAME: ['bin/scistree'],  # Corrected path
    },
    include_package_data=True,
)
