from setuptools import setup, find_packages

setup(
    name='magus',
    version='0.1.0',
    packages=find_packages(),  # This will include 'magus' and all sub-packages
    entry_points={
        'console_scripts': [
            'magus=magus.magus_main:main',  # Ensure the entry point specifies the package path
        ],
    },
    include_package_data=True,
    install_requires=[
        # List any required dependencies for Python here, or leave it empty if using Conda exclusively
    ],
)

