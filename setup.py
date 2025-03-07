from setuptools import setup, find_packages

setup(
    name='hic_basic',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        # 列出所有依赖
    ],
    entry_points={
        'console_scripts': [
            'hic_basic=hic_basic.__main__:main',
        ],
    },
)