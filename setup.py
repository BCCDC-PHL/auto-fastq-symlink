from setuptools import setup, find_namespace_packages


setup(
    name='auto-fastq-symlink',
    version='0.1.0-alpha',
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "auto-fastq-symlink = auto_fastq_symlink.__main__:main",
        ]
    },
    scripts=[],
    package_data={
        "auto_fastq_symlink.resources": ["*.json"],
    },
    install_requires=[
        "sqlalchemy==1.4.40",
        "alembic==1.8.0",
        "jsonschema",
        "uvicorn",
        "fastapi",
    ],
    description=' Automated symlinking of sequence data',
    url='https://github.com/BCCDC-PHL/auto-fastq-symlink',
    author='Dan Fornika',
    author_email='dan.fornika@bccdc.ca',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
