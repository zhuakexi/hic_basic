from setuptools import setup, find_packages

setup(
    name='hic_basic',
    version='0.1',
    packages=find_packages(),
    package_data={
        'hic_basic': ['data/*.json', 'data/viewer/*'],
    },
    install_requires=[
        # Core dependencies
        'numpy',
        'pandas',
        'scipy',
        'matplotlib',
        'plotly',
        'kaleido',
        'cooler',
        'cooltools',
        'scikit-learn',
        'seaborn',
        # BigWig support
        'pyBigWig',
    ],
    entry_points={
        'console_scripts': [
            'hic_basic=hic_basic.__main__:main',
        ],
    },
)
