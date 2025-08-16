import os
from setuptools import setup

def get_version_number():
  path_to_written_version = 'eccp/__init__.py'
  with open(path_to_written_version) as initPY:
    for line in initPY:
      if line.startswith('__version__'):
        version = eval(line.rstrip().replace('__version__ = ',''))
        break
  return version

def get_long_description():
  this_directory = os.path.abspath(os.path.dirname(__file__))
  readme_file = os.path.join(this_directory, 'README.md')
  with open(readme_file, encoding='utf-8') as f:
    long_description = f.read()
  return long_description

def find_packages(root):
  packages = []
  for root, dirs, files in os.walk(root, topdown=False):
    for name in files:
      if name.endswith('.py'):
        file = os.path.relpath(os.path.join(root, name))
        dirname = os.path.dirname(file)
        rel_dirname = os.path.relpath(dirname)
        if not rel_dirname in packages:
          packages.append(rel_dirname)
  return sorted(packages)

def find_scripts():
  scripts = []
  scripts_folders = []
  for scripts_folder in scripts_folders:
    for root, dirs, files in os.walk('eccp/'+scripts_folder, topdown=False):
      for file in files:
        if file.endswith('.py') and not 'Main' in file:
          filepath = os.path.relpath(os.path.join(root, file))
          scripts.append(filepath)
      dirs[:] = []
  return sorted(scripts)

setup(name='ECCP',
      packages=find_packages(root='eccp'),
      scripts=['bin/eccp']+find_scripts(),
      version=get_version_number(),
      description="This program takes a crystal structure (of a OPV or any crystal) and will separate the crystal into its individual molecules, as well as provide the dimers of these molecules. This program will then provide Gaussian input files for performing atomic transition charges upon the unique molecules, as well as for performing excitation energy transfer calculations upon the unique dimers.",
      long_description=get_long_description(),
      long_description_content_type='text/markdown',
      author='Dr. Geoffrey R. Weal and Dr. Anna L. Garden',
      author_email='anna.garden@otago.ac.nz',
      url = 'https://blogs.otago.ac.nz/annagarden/',
      download_url = 'https://github.com/GardenGroupUO/Electronic_Crystal_Calculation_Prep/archive/v'+str(get_version_number())+'.tar.gz',
      license='GNU AFFERO GENERAL PUBLIC LICENSE',
      zip_safe=False,
      keywords = ['victoria-university', 'victoria-university-of-wellington', 'university-of-wellington', 'wellington-university', 'atomic-simulation-environment', 'organic-photovoltaics', 'OPV'],
      install_requires=['numpy', 'ase>=3.19.0', 'packaging', 'networkx', 'tqdm', 'xlsxwriter'],
      classifiers=[
        'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Science/Research',      # Define that your audience are developers
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Natural Language :: English',
        'License :: OSI Approved :: GNU Affero General Public License v3',   # Again, pick a license
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        ],
      )

