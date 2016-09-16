import os
from fluxpart.fluxpart import flux_partition

basedir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
with open(os.path.join(basedir, 'VERSION')) as vf:
    __version__ = vf.read().strip()
