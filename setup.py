from setuptools import setup, find_packages
setup(
    name='chemi',
    version='0.1.0',
    install_requires=["rdkit"],
    packages=find_packages(include=['chemi', 'chemi.*'])
)