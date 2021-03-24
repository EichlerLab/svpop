__all__ = []

from . import anno
from . import plot
from . import tracks

from . import callerset
from . import gt
from . import pd
from . import ref
from . import refseq
from . import rules
from . import sampleset
from . import seq
from . import sm
from . import svlenoverlap
from . import svmerge
from . import svset
from . import util
from . import varbed
from . import variant
from . import variant_extern
from . import varset
from . import vcf

#import pkgutil
#
#__all__ = []
#
## Using method described in:
## https://stackoverflow.com/questions/3365740/how-to-import-all-submodules
#
# for loader, mod_name, is_pkg in pkgutil.walk_packages(__path__):
#     __all__.append(mod_name)
#     _module = loader.find_module(mod_name).load_module(mod_name)
#     globals()[mod_name] = _module
