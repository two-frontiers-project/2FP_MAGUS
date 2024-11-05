from setuptools import setup, find_packages

setup(
    name='magus',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'magus=magus.magus_main:main',
        ],
    },
    data_files=[
        ('bin', [
            'bin/magus_assembly_helper.sh', 
            'bin/sorenson-g',
            'bin/shi7_trimmer',
            'bin/minigzip',
            'bin/akmer100b',
            'bin/lingenome',
            'bin/fac',
            'bin/bestmag2',
            'bin/spamw2',
            'bin/megahit-g',
            'bin/canolax5',
        ]),
    ],
)

