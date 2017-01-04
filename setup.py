try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'My Project',
    'author': 'Karl Lapo',
    'url': 'github.com/klapo',
    'author_email': 'lapo.karl@gmail.com',
    'version': '0.1',
#    'install_requires': ['nose'],
    'packages': ['pyir'],
    'scripts': [],
    'name': 'pyir'
}

setup(**config)
