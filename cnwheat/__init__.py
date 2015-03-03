# -*- coding: latin-1 -*-
"""
    cnwheat
    ~~~~~~~

    The model CN-Wheat.
    
    CN-Wheat computes the CN exchanges in a wheat architecture. See:
    
        * :mod:`cnwheat.simulation` for the front-end of the model,
        * :mod:`cnwheat.model` for the equations of the model,
        * :mod:`cnwheat.parameters` for the parameters of the model.
         
    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
    
    .. seealso:: Barillot et al. 2014.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

__version__  = '0.0.1'

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())