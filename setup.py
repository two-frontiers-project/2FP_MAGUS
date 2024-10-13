from setuptools import setup, find_packages

setup(
    name='magus',
    version='0.1.0',
    packages=find_packages(),  # Automatically find and include all sub-packages
    install_requires=[
        # Specify any Python dependencies here (you can leave this empty if handled by Conda)
    ],
    entry_points={
        'console_scripts': [
            'magus=magus_main:main',  # Register the CLI command 'magus'
        ],
    },
    include_package_data=True,
)

