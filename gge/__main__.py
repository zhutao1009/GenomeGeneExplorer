"""`python -m gge` 入口。"""

import sys

from .ui.app import main

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
