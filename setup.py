from setuptools import setup

setup(
    name='GMMchi_scs',
    version='0.1',
    author='Jeff Liu',
    author_email='jeffliu6068@gmail.com',
    packages=['GMMchi_scs'],
    url='https://github.com/jeffliu6068/GMMchi_scs',
    license='LICENSE.txt',
    description='Pre-print version of the GMMchi-based single cell postprocessing pipeline',
    long_description=open('README.md').read(),
    install_requires=[
        'pandas',
        'scipy',
        'numpy',
        'matplotlib',
        'seaborn',
        'sklearn',
        'tqdm',
        'umap-learn',
        'scprep'
    ],
)
