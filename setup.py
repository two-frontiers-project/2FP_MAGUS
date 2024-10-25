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
    install_requires=[
        # List any required dependencies
    ],
)

