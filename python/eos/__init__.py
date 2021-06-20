# Copyright (c) 2018 Frederik Beaujean
# Copyright (c) 2017, 2018, 2020 Danny van Dyk
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

from .config import *

# make sure that EOS_HOME points to the location of the wheel supplied data
# if unset.
import os as _os
try:
    if is_wheel:
        if not 'EOS_HOME' in _os.environ:
            _os.environ['EOS_HOME'] = _os.path.normpath(_os.path.join(_os.path.dirname(__file__), '..', '_eos_data/'))
except NameError:
    pass

from _eos import *
from _eos import _register_callback
from .data import *
from .plot import *
from .analysis import Analysis, BestFitPoint
from .analysis_file import AnalysisFile
from .constraint import Constraints
from .observable import Observables
from .parameter import Parameters
from .reference import References
from .signal_pdf import SignalPDF, SignalPDFs

import logging
logger = logging.getLogger('EOS')
logger.setLevel(logging.INFO)

def debug(msg, *args, **kwargs):
    logger.debug(msg, *args, **kwargs)

def error(msg, *args, **kwargs):
    logger.error(msg, *args, **kwargs)

def info(msg, *args, **kwargs):
    logger.info(msg, *args, **kwargs)

def warn(msg, *args, **kwargs):
    logger.warn(msg, *args, **kwargs)

def _log_callback(id, level, msg):
    logger.info("{id} {level} {msg}".format(id=id, level=level, msg=msg))

_register_callback(_log_callback)

import time as _time
import os as _os
def installation_time():
    return _time.ctime(_os.path.getmtime(eos.__file__))

def installation_dir():
    return _os.path.dirname(eos.__file__)

__ipython__ = False
try:
    if __IPYTHON__:
        __ipython__ = True
        ip = get_ipython()
        html_formatter = ip.display_formatter.formatters['text/html']

        from .ipython import __format_Parameter, __format_KinematicVariable, __format_Kinematics, __format_Options, __format_Observable, __format_ObservableEntry, __format_GoodnessOfFit
        html_formatter.for_type(Parameter, __format_Parameter)
        html_formatter.for_type(KinematicVariable, __format_KinematicVariable)
        html_formatter.for_type(Kinematics, __format_Kinematics)
        html_formatter.for_type(Options, __format_Options)
        html_formatter.for_type(Observable, __format_Observable)
        html_formatter.for_type(ObservableEntry, __format_ObservableEntry)
        html_formatter.for_type(GoodnessOfFit, __format_GoodnessOfFit)

        # print logger in notebook
        import logging
        logger = logging.getLogger('EOS')
        logger.setLevel(logging.INFO)
        logging.basicConfig(stream=sys.stderr, level=logging.INFO)
except NameError as e:
    pass
