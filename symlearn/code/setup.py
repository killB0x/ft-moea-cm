from setuptools import setup

setup(
    name='ft_learning',
    version='0.2',
    author='L.A. Jimenez-Roa, M. Volk',
    author_email='m.volk@utwente.nl',
    maintainer='M. Volk',
    maintainer_email='m.volk@utwente.nl',
    description='ft_learning - Learning Fault Trees from Data',
    zip_safe=False,
    install_requires=['pandas',  # Handling data sets
                      'numpy',  # Handling data sets
                      'sympy',  # Simplification of Boolean formulas
                      'pyeda',  # Simplification of Boolean formulas
                      'scipy',  # Export/Import to matlab
                      'sortedcontainers',  # Sorted set
                      ],
    python_requires='>=3',
)
