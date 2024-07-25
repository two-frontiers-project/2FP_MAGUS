from setuptools import setup, find_packages

setup(
    name='magus',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
    ],
    entry_points={
        'console_scripts': [
            'magus-taxonomy=magus.taxonomy:main',
            'magus-assemble=magus.assemble:main',
            'magus-resolvegenomes=magus.resolvegenomes:main',
            'magus-genelevel=magus.genelevel:main',
        ],
    },
)

