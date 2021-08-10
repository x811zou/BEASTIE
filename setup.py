from setuptools import setup

with open("requirements.txt", "r") as requirements_file:
	install_requirements = requirements_file.readlines()

setup(
	install_requires = install_requirements,
)
