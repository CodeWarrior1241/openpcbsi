"""The board generator lives in the package (openpcbsi/boardgen.py); the tests import it from here."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from openpcbsi.boardgen import *  # noqa: E402,F401,F403
from openpcbsi.boardgen import Board, LAYER, STACKUPS  # noqa: E402,F401
