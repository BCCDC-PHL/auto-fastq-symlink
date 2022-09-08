from setuptools import setup, find_packages


setup(
    name='auto-fastq-symlink',
    version='0.1.0',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "auto-fastq-symlink = auto_fastq_symlink.__main__:main",
        ]
    },
    scripts=[],
    package_data={},
    install_requires=[
        "sqlalchemy==1.4.40",
        "alembic==1.8.0",
        "uvicorn",
        "fastapi",
    ],
    description='',
    url='',
    author='',
    author_email='',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
