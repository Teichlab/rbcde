from setuptools import setup, find_packages

setup(
	name='rbcde',
	version='1.0.0',
	description='Rank-biserial correlation',
	url='https://github.com/Teichlab/rbcde',
	packages=find_packages(exclude=['docs', 'examples']),
	install_requires=['numpy','scipy','pandas'],
	author='Krzysztof Polanski',
	author_email='kp9@sanger.ac.uk',
	license='MIT'
)