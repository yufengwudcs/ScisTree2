from setuptools import setup, find_packages, Distribution
from setuptools.command.build_py import build_py
import subprocess
import shutil
import os
import sys
from pathlib import Path


PACKAGE_NAME = 'scistree2'
# PROJECT_ROOT = str(Path(__file__).parent.resolve())
PACKAGE_DIR = 'scistree2'
SRC_DIR = "src"
CPP_DIR = "src"
# CPP_DIR = f'{PROJECT_ROOT}/src'
# PACKAGE_DIR = f'{PROJECT_ROOT}/scistree2'


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True  # <--- The magic line. Forces "platlib" (platform-specific) wheel.
    

# class BuildWithMake(build_py):
#     def run(self):
#         self.build_cpp_with_make()
#         super().run()
    
#     def build_cpp_with_make(self):
#         print("Building C++ project with make...")
#         try:
#             subprocess.run(['make'], cwd=CPP_DIR, check=True)
#         except subprocess.CalledProcessError as e:
#             print("Make failed.")
#             raise e
        
#         bin_dir = os.path.join(PACKAGE_DIR, 'bin')
#         os.makedirs(bin_dir, exist_ok=True)
#         shutil.copy(
#             os.path.join(CPP_DIR, 'scistree'),
#             os.path.join(bin_dir, 'scistree')
#         )
class CustomBuild(build_py):
    def run(self):
        if os.environ.get("CIBUILDWHEEL") == "1":
            print("--- Detected CIBUILDWHEEL: Skipping setup.py compilation (relying on pre-built binary) ---")
        else:
            # Local install (pip install .): We must compile.
            self.compile_cpp()
        super().run()

    def compile_cpp(self):
        print(f"--- Compiling C++ extensions for {sys.platform} ---")
        
        bin_dir = os.path.join(PACKAGE_DIR, "bin")
        os.makedirs(bin_dir, exist_ok=True)

        if sys.platform == "win32":
            # Windows: Use nmake
            build_cmd = "nmake /f Makefile.win"
            clean_cmd = "nmake /f Makefile.win clean"
            binary_name_src = "scistree.exe"
            binary_name_dest = "scistree.exe"
        else:
            # Linux/Mac: Use standard make
            build_cmd = "make"
            clean_cmd = "make clean"
            binary_name_src = "scistree"
            binary_name_dest = "scistree"

        try:
            # Run clean then build inside the 'src' directory
            subprocess.check_call(clean_cmd, cwd=SRC_DIR, shell=True)
            subprocess.check_call(build_cmd, cwd=SRC_DIR, shell=True)
        except subprocess.CalledProcessError as e:
            print("Error during compilation:", e)
            sys.exit(1)
        except FileNotFoundError:
            print(f"Error: Build tool not found. Ensure '{build_cmd[0]}' is in PATH.")
            sys.exit(1)

        src_path = os.path.join(SRC_DIR, binary_name_src)
        dest_path = os.path.join(bin_dir, binary_name_dest)
        
        if os.path.exists(src_path):
            print(f"Moving {src_path} -> {dest_path}")
            shutil.copy(src_path, dest_path) 
        else:
            print(f"Error: Compiled binary {binary_name_src} not found in {SRC_DIR}")
            sys.exit(1)

setup(
    name=PACKAGE_NAME,
    packages=[PACKAGE_DIR],
    cmdclass={'build_py': CustomBuild},
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
