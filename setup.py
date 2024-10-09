from setuptools import setup, find_packages

setup(
    name='magus',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        # Add dependencies here, e.g., 'numpy', 'pandas'
    ],
    entry_points={
        'console_scripts': [
            'magus=magus.qc:main',
            'magus=magus.single_assembly:main',
            'magus=magus.single_binning:main',
            'magus=magus.cluster_contigs:main',
            'magus=magus.coassembly_binning:main',
            'magus=magus.coassembly:main',
            'magus=magus.taxonomy:main',
        ],
    },
)

