import logging

from main import __version__

logger = logging.getLogger(__name__) 
logger.addHandler(logging.NullHandler())
